#!/usr/bin/env python3
"""Build the SQLite database behind the Phallusia expression app.

Inputs (all produced by the pipeline, see ../README.md):
  --gtf       annotation_with_MT.gtf.gz        RefSeq GTF + mitochondrial genes
  --genome    genome_with_MT.fna.gz            RefSeq genome + mitochondrion
  --rna       transcripts.fna.gz               RefSeq mRNA/ncRNA sequences     (rna.fna.gz)
  --cds       cds.fna.gz                       RefSeq CDS sequences            (cds_from_genomic.fna.gz)
  --protein   proteins.faa.gz                  RefSeq proteins
  --tpm       salmon.merged.gene_tpm.tsv       merged gene x sample TPM (03d_merge_quant.sbatch; same layout as nf-core/rnaseq)
  --counts    salmon.merged.gene_counts.tsv    merged gene x sample counts (same layout)
  --samples   sample_metadata.csv              sample, run_accession, ena_sample, stage, hotta_stage, replicate, dev_order
  --eggnog    phmamm.emapper.annotations       eggNOG-mapper output (optional)
  --extra-protein MT_OZ281502.1.proteins.faa  more proteins, e.g. the MITOS2 mitochondrial ones (optional, repeatable;
                                               ids = protein_id of the GTF CDS lines; 01c_annotate_mt_mitos.sh)
  --mt-prefix MT-                              genes with this prefix (the mitochondrial genome) get their own TPM scale.
                                               Per sample, Salmon's TPM is rescaled so that the nuclear genes sum to 10^6
                                               and, separately, the mitochondrial genes WITHOUT the 2 rRNAs sum to 10^6
                                               (16S rRNA alone would be ~99% of the MT total). The MT rRNAs are scaled
                                               with the same factor, so they can exceed 10^6. Why: MT genes are
                                               transcribed in the mitochondria, not the nucleus, and take 15-30% of all TPM
                                               (changing with stage), which would shift every nuclear gene between stages.
                                               --mt-prefix "" stores Salmon's TPM unchanged. Counts are never changed.

Output: one SQLite file (--out). Standard library only; needs ~2 GB RAM for the real genome.
"""
import argparse
import csv
import datetime
import gzip
import re
import sqlite3
import sys
import zlib

ATTR = re.compile(r'(\w+) "([^"]*)"')
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
VERSION = "0.2"


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, buf = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)


def parse_gtf(path):
    genes, tx = {}, {}
    for line in opener(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        seqid, _, feat, s, e, _, strand, _, attr = f[:9]
        pairs = ATTR.findall(attr)
        a = dict(pairs)
        gid = a.get("gene_id")
        if not gid:
            continue
        s, e = int(s), int(e)
        if feat == "gene":
            g = genes.setdefault(gid, {})
            geneid = next((v.split(":", 1)[1] for k, v in pairs if k == "db_xref" and v.startswith("GeneID:")), "")
            g.update(seqid=seqid, start=s, end=e, strand=strand, biotype=a.get("gene_biotype", ""),
                     description=a.get("description", ""), ncbi_gene_id=geneid,
                     symbol=a.get("gene_name") or a.get("gene") or gid)
            continue
        tid = a.get("transcript_id")
        if not tid:
            continue
        t = tx.setdefault(tid, {"gene_id": gid, "seqid": seqid, "strand": strand, "product": "", "biotype": "",
                                "exons": [], "cds": [], "protein_id": ""})
        if feat == "transcript":
            t["product"] = a.get("product", t["product"])
            t["biotype"] = a.get("transcript_biotype", a.get("gene_biotype", t["biotype"]))
        elif feat == "exon":
            t["exons"].append((s, e))
            t["product"] = t["product"] or a.get("product", "")
        elif feat == "CDS":
            t["cds"].append((s, e))
            t["protein_id"] = a.get("protein_id", t["protein_id"])
            t["product"] = t["product"] or a.get("product", "")
    for tid, t in tx.items():                       # genes without a 'gene' line (not expected for RefSeq)
        g = genes.setdefault(t["gene_id"], {})
        if "seqid" not in g:
            g.update(seqid=t["seqid"], strand=t["strand"], biotype=t["biotype"], description=t["product"],
                     ncbi_gene_id="", symbol=t["gene_id"], start=min(x[0] for x in t["exons"]), end=max(x[1] for x in t["exons"]))
    return genes, tx


def join_segments(genome, seqid, strand, segs):
    seq = "".join(genome[seqid][s - 1:e] for s, e in sorted(segs))
    return revcomp(seq) if strand == "-" else seq


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("gtf", "genome", "rna", "cds", "protein", "tpm", "counts", "samples", "out"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--eggnog")
    ap.add_argument("--extra-protein", action="append", default=[])
    ap.add_argument("--mt-prefix", default="MT-")
    ap.add_argument("--demo", action="store_true", help="mark the database as simulated demo data")
    ap.add_argument("--assembly", default="GCF_965637545.1")
    ap.add_argument("--annotation", default="RefSeq annotation RS_2026_02")
    ap.add_argument("--reads", default="ENA PRJEB25149 (E-MTAB-6528)")
    ap.add_argument("--pipeline", default="STAR (two-pass) + Salmon (alignment-based mode)")
    a = ap.parse_args()
    warn = []

    print("Reading annotation ...")
    genes, tx = parse_gtf(a.gtf)
    print(f"  {len(genes)} genes, {len(tx)} transcripts")

    print("Reading genome ...")
    genome = {n.split()[0]: s.upper() for n, s in read_fasta(a.genome)}
    missing = {g["seqid"] for g in genes.values()} - set(genome)
    if missing:
        sys.exit(f"ERROR: GTF sequences missing from the genome FASTA: {sorted(missing)[:5]}")

    print("Reading transcript / CDS / protein sequences ...")
    rna = {n.split()[0]: s.upper() for n, s in read_fasta(a.rna)}
    cds = {}
    for n, s in read_fasta(a.cds):
        m = re.search(r"\[protein_id=([^\]]+)\]", n)
        if m:
            cds[m.group(1)] = s.upper()
    prot = {n.split()[0]: s for n, s in read_fasta(a.protein)}
    for path in a.extra_protein:
        extra = {n.split()[0]: s for n, s in read_fasta(path)}
        clash = set(extra) & set(prot)
        if clash:
            sys.exit(f"ERROR: protein ids of {path} also in {a.protein}: {sorted(clash)[:3]}")
        prot.update(extra)
        print(f"  {len(extra)} extra proteins from {path}")

    # ---- samples ---------------------------------------------------------------------------
    samples = sorted(csv.DictReader(open(a.samples)), key=lambda r: (int(r["dev_order"]), int(r["replicate"])))
    for i, r in enumerate(samples):
        r["sample_idx"] = i
        r["stage_label"] = f'{r["stage"]} ({r["hotta_stage"]})' if r.get("hotta_stage") else r["stage"]
    sample_idx = {r["sample"]: r["sample_idx"] for r in samples}

    db = sqlite3.connect(a.out)
    db.executescript("""
        PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;
        DROP TABLE IF EXISTS meta; DROP TABLE IF EXISTS samples; DROP TABLE IF EXISTS genes; DROP TABLE IF EXISTS transcripts;
        DROP TABLE IF EXISTS seq_tx; DROP TABLE IF EXISTS seq_gene; DROP TABLE IF EXISTS expression; DROP TABLE IF EXISTS eggnog;
        DROP TABLE IF EXISTS search; DROP TABLE IF EXISTS identifiers;
        CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
        CREATE TABLE samples(sample_idx INTEGER PRIMARY KEY, sample TEXT, run_accession TEXT, ena_sample TEXT, stage TEXT,
                             hotta_stage TEXT, stage_label TEXT, replicate INTEGER, dev_order INTEGER);
        CREATE TABLE genes(gene_id TEXT PRIMARY KEY, symbol TEXT, description TEXT, biotype TEXT, seqid TEXT, start INTEGER,
                           end INTEGER, strand TEXT, ncbi_gene_id TEXT, n_transcripts INTEGER, primary_transcript TEXT,
                           peak_stage TEXT, peak_dev_order INTEGER, peak_tpm REAL);
        CREATE TABLE transcripts(transcript_id TEXT PRIMARY KEY, gene_id TEXT, product TEXT, biotype TEXT, protein_id TEXT,
                                 n_exons INTEGER, mrna_len INTEGER, cds_len INTEGER, aa_len INTEGER, is_primary INTEGER);
        CREATE TABLE seq_tx(transcript_id TEXT PRIMARY KEY, mrna_z BLOB, cds_z BLOB, protein TEXT);  -- DNA zlib-compressed
        CREATE TABLE seq_gene(gene_id TEXT PRIMARY KEY, length INTEGER, seq_zlib BLOB);
        CREATE TABLE expression(gene_id TEXT, sample_idx INTEGER, tpm REAL, counts REAL, PRIMARY KEY(gene_id, sample_idx)) WITHOUT ROWID;
        CREATE TABLE eggnog(gene_id TEXT PRIMARY KEY, protein_id TEXT, description TEXT, preferred_name TEXT, cog_category TEXT,
                            og TEXT, max_annot_lvl TEXT, gos TEXT, ec TEXT, kegg_ko TEXT, kegg_pathway TEXT, pfams TEXT,
                            seed_ortholog TEXT, evalue TEXT, score TEXT);
        CREATE TABLE search(gene_id TEXT PRIMARY KEY, symbol_l TEXT, desc_l TEXT, other_l TEXT);
        CREATE TABLE identifiers(identifier_l TEXT, gene_id TEXT);
    """)
    db.executemany("INSERT INTO samples VALUES(?,?,?,?,?,?,?,?,?)",
                   [(r["sample_idx"], r["sample"], r["run_accession"], r.get("ena_sample", ""), r["stage"], r.get("hotta_stage", ""),
                     r["stage_label"], int(r["replicate"]), int(r["dev_order"])) for r in samples])

    # ---- transcripts and sequences --------------------------------------------------------------
    print("Building transcript and gene sequences ...")
    n_rebuilt, n_odd = 0, 0       # n_odd: rebuilt although RefSeq's RNA FASTA should have them (not tRNA, not MT)
    by_gene = {}
    seq_rows, tx_rows = [], []
    for tid, t in tx.items():
        mrna = rna.get(tid)
        if mrna is None:
            mrna = join_segments(genome, t["seqid"], t["strand"], t["exons"])
            n_rebuilt += 1
            n_odd += not (t["biotype"] == "tRNA" or t["gene_id"].startswith("MT-"))
        c = cds.get(t["protein_id"], "") if t["protein_id"] else ""
        if t["cds"] and not c:
            c = join_segments(genome, t["seqid"], t["strand"], t["cds"])
            if not t["gene_id"].startswith("MT-"):   # MT CDS are never in the RefSeq CDS FASTA
                warn.append(f"CDS of {tid} rebuilt from the genome (not in the CDS FASTA)")
        p = prot.get(t["protein_id"], "") if t["protein_id"] else ""
        t.update(mrna_len=len(mrna), cds_len=len(c), aa_len=len(p))
        by_gene.setdefault(t["gene_id"], []).append(tid)
        seq_rows.append((tid, zlib.compress(mrna.encode(), 6), zlib.compress(c.encode(), 6) if c else None, p))
    print(f"  {n_rebuilt} transcripts were not in the RefSeq RNA FASTA; their sequence was rebuilt from exons "
          f"({n_rebuilt - n_odd} tRNAs / mitochondrial genes, which RefSeq does not put there; {n_odd} others)")
    if n_odd > 0.05 * max(len(tx), 1):
        warn.append(f"{n_odd} of {len(tx)} transcripts (not tRNA/MT) had to be rebuilt from exons; is the RNA FASTA from the same release?")
    for gid, tids in by_gene.items():
        best = max(tids, key=lambda x: (tx[x]["cds_len"], tx[x]["mrna_len"]))
        genes[gid]["primary_transcript"] = best
        for tid in tids:
            tx[tid]["is_primary"] = int(tid == best)
    for tid, t in tx.items():
        tx_rows.append((tid, t["gene_id"], t["product"], t["biotype"], t["protein_id"], len(t["exons"]),
                        t["mrna_len"], t["cds_len"], t["aa_len"], t.get("is_primary", 0)))
    db.executemany("INSERT INTO transcripts VALUES(?,?,?,?,?,?,?,?,?,?)", tx_rows)
    db.executemany("INSERT INTO seq_tx VALUES(?,?,?,?)", seq_rows)
    del seq_rows

    gene_seq_rows = []
    for gid, g in genes.items():
        s = genome[g["seqid"]][g["start"] - 1:g["end"]]
        if g["strand"] == "-":
            s = revcomp(s)
        gene_seq_rows.append((gid, len(s), zlib.compress(s.encode(), 6)))
    db.executemany("INSERT INTO seq_gene VALUES(?,?,?)", gene_seq_rows)
    del gene_seq_rows, genome

    # ---- expression ----------------------------------------------------------------------------
    print("Reading expression tables ...")

    def read_matrix(path):
        with open(path) as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            cols = hdr[2:]
            unknown = [c for c in cols if c not in sample_idx]
            if unknown:
                sys.exit(f"ERROR: columns of {path} not in sample_metadata.csv: {unknown[:5]}")
            return cols, {row[0]: [float(x) if x not in ("", "NA") else 0.0 for x in row[2:]] for row in rd}

    cols, tpm = read_matrix(a.tpm)
    cols2, cnt = read_matrix(a.counts)
    if cols != cols2:
        sys.exit("ERROR: TPM and counts tables have different sample columns")
    not_in_gtf = [g for g in tpm if g not in genes]
    if len(not_in_gtf) > 0.01 * len(tpm):
        sys.exit(f"ERROR: {len(not_in_gtf)} genes of the expression table are not in the GTF (e.g. {not_in_gtf[:3]}). Same annotation used for the run?")
    if not_in_gtf:
        warn.append(f"{len(not_in_gtf)} expression rows skipped: gene not in GTF ({not_in_gtf[:3]} ...)")
    tpm_note = "Salmon TPM"
    mt_ids = [g for g in tpm if a.mt_prefix and g.startswith(a.mt_prefix)]
    if mt_ids:
        mt_set = set(mt_ids)
        mt_rrna = {g for g in mt_ids if g in genes and genes[g].get("biotype") == "rRNA"}
        all_mt_sum = [sum(tpm[g][j] for g in mt_ids) for j in range(len(cols))]
        mt_sum = [sum(tpm[g][j] for g in mt_ids if g not in mt_rrna) for j in range(len(cols))]
        nuc_sum = [sum(v[j] for g, v in tpm.items() if g not in mt_set) for j in range(len(cols))]
        if min(mt_sum) <= 0 or min(nuc_sum) <= 0:
            sys.exit("ERROR: a sample has no TPM in the mitochondrial or in the nuclear genes; check --mt-prefix")
        for gid in tpm:
            tot = mt_sum if gid in mt_set else nuc_sum
            tpm[gid] = [v * 1e6 / t for v, t in zip(tpm[gid], tot)]
        tpm_note = (f"nuclear and mitochondrial ({a.mt_prefix}...) genes are on separate TPM scales: per sample, Salmon's TPM is "
                    f"rescaled so that the nuclear genes sum to 10^6 and the mitochondrial genes except the rRNAs "
                    f"({', '.join(sorted(mt_rrna)) or 'none'}) sum to 10^6; the MT rRNAs use the same factor and can exceed 10^6")
        print(f"  TPM: {len(mt_ids)} mitochondrial genes put on their own scale (without the rRNAs {', '.join(sorted(mt_rrna)) or 'none'}); "
              f"all MT genes took {min(all_mt_sum)/1e4:.1f}-{max(all_mt_sum)/1e4:.1f}% of Salmon's TPM, the non-rRNA ones "
              f"{min(mt_sum)/1e4:.2f}-{max(mt_sum)/1e4:.2f}%")
    rows = []
    dev_of = {r["sample_idx"]: int(r["dev_order"]) for r in samples}
    peak = {}
    for gid, vals in tpm.items():
        if gid not in genes:
            continue
        per_stage = {}
        for c, v, n in zip(cols, vals, cnt[gid]):
            si = sample_idx[c]
            rows.append((gid, si, v, n))
            per_stage.setdefault(dev_of[si], []).append(v)
        means = {k: sum(x) / len(x) for k, x in per_stage.items()}
        k = max(means, key=means.get)
        peak[gid] = (k, means[k])
    db.executemany("INSERT INTO expression VALUES(?,?,?,?)", rows)
    label_of = {int(r["dev_order"]): r["stage_label"] for r in samples}
    no_expr = [g for g in genes if g not in peak]
    if no_expr:
        warn.append(f"{len(no_expr)} genes of the GTF have no expression row")
    print(f"  {len(rows)} expression values for {len(peak)} genes")

    # ---- eggNOG --------------------------------------------------------------------------------
    egg_by_gene = {}
    prot2gene = {t["protein_id"]: t["gene_id"] for t in tx.values() if t["protein_id"]}
    if a.eggnog:
        print("Reading eggNOG-mapper annotations ...")
        names, n_q = None, 0
        with open(a.eggnog) as fh:
            for line in fh:
                if line.startswith("#query"):
                    names = line.lstrip("#").rstrip("\n").split("\t")
                    continue
                if line.startswith("#") or not line.strip() or names is None:
                    continue
                r = dict(zip(names, line.rstrip("\n").split("\t")))
                gid = prot2gene.get(r["query"])
                if gid is None:
                    continue
                n_q += 1
                egg_by_gene.setdefault(gid, []).append(r)
        if names is None:
            sys.exit("ERROR: no '#query' header line in the eggNOG file")
        if n_q == 0:
            warn.append("no eggNOG query matched a RefSeq protein id; was eggNOG run on proteins.faa from this release?")
        clean = lambda x: "" if x in (None, "-") else x
        egg_rows = []
        for gid, rs in egg_by_gene.items():
            prim = tx[genes[gid]["primary_transcript"]]["protein_id"]
            r = next((x for x in rs if x["query"] == prim), None) or max(rs, key=lambda x: float(clean(x.get("score")) or 0))
            egg_rows.append((gid, r["query"], clean(r.get("Description")), clean(r.get("Preferred_name")), clean(r.get("COG_category")),
                             clean(r.get("eggNOG_OGs")), clean(r.get("max_annot_lvl")), clean(r.get("GOs")), clean(r.get("EC")),
                             clean(r.get("KEGG_ko")), clean(r.get("KEGG_Pathway")), clean(r.get("PFAMs")),
                             clean(r.get("seed_ortholog")), clean(r.get("evalue")), clean(r.get("score"))))
        db.executemany("INSERT INTO eggnog VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", egg_rows)
        print(f"  eggNOG annotation for {len(egg_rows)} genes")
        egg = {r[0]: r for r in egg_rows}
    else:
        egg = {}

    # ---- genes, search, identifiers ----------------------------------------------------------------
    gene_rows, search_rows, id_rows = [], [], []
    for gid, g in genes.items():
        tids = by_gene.get(gid, [])
        desc = g.get("description") or next((tx[t]["product"] for t in tids if tx[t]["product"]), "")
        if not desc and gid.startswith("MT-"):
            desc = f"mitochondrial gene {g['symbol']}"
        pk = peak.get(gid)
        gene_rows.append((gid, g["symbol"], desc, g["biotype"], g["seqid"], g["start"], g["end"], g["strand"], g.get("ncbi_gene_id", ""),
                          len(tids), g.get("primary_transcript", ""), label_of.get(pk[0]) if pk else None, pk[0] if pk else None,
                          pk[1] if pk else None))
        e = egg.get(gid)
        prods = " ".join(sorted({tx[t]["product"] for t in tids if tx[t]["product"]}))
        other = " ".join(x for x in [prods, e[2] if e else "", e[3] if e else "", e[11].replace(",", " ") if e else "",
                                     e[7].replace(",", " ") if e else "", e[9] if e else "", g.get("ncbi_gene_id", ""),
                                     " ".join(tids), " ".join(tx[t]["protein_id"] for t in tids if tx[t]["protein_id"])] if x)
        search_rows.append((gid, f'{gid} {g["symbol"]}'.lower(), desc.lower(), other.lower()))
        ids = {gid, g["symbol"], g.get("ncbi_gene_id", ""), "geneid:" + g.get("ncbi_gene_id", "")} - {"", "geneid:"}
        for t in tids:
            ids.add(t)
            ids.add(t.split(".")[0])
            pid = tx[t]["protein_id"]
            if pid:
                ids.add(pid)
                ids.add(pid.split(".")[0])
        id_rows += [(i.lower(), gid) for i in ids]
    db.executemany("INSERT INTO genes VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)", gene_rows)
    db.executemany("INSERT INTO search VALUES(?,?,?,?)", search_rows)
    db.executemany("INSERT INTO identifiers VALUES(?,?)", id_rows)
    db.execute("CREATE INDEX idx_tx_gene ON transcripts(gene_id)")
    db.execute("CREATE INDEX idx_ident ON identifiers(identifier_l)")

    meta = {"title": "Phallusia mammillata expression atlas", "assembly": a.assembly, "annotation": a.annotation,
            "reads": a.reads, "pipeline": a.pipeline, "eggnog": "eggNOG-mapper v2 on RefSeq proteins" if a.eggnog else "none",
            "built": datetime.date.today().isoformat(), "builder_version": VERSION, "is_demo": "1" if a.demo else "0",
            "n_genes": str(len(genes)), "n_samples": str(len(samples)), "tpm_note": tpm_note}
    db.executemany("INSERT INTO meta VALUES(?,?)", meta.items())
    db.commit()
    db.execute("VACUUM")
    db.close()
    print(f"Wrote {a.out}")
    for w in warn[:20]:
        print("WARNING:", w)
    if len(warn) > 20:
        print(f"... and {len(warn) - 20} more warnings")


if __name__ == "__main__":
    main()
