#!/usr/bin/env python3
"""Make SIMULATED input files in the same formats as the real ones, so the database builder and
the Shiny app can be tested end to end without the real data:

    genome.fna.gz  annotation.gtf.gz  rna.fna.gz  cds.fna.gz  proteins.faa.gz   (RefSeq-style)
    salmon.merged.gene_tpm.tsv  salmon.merged.gene_counts.tsv                   (merged-table layout of the pipeline)
    phmamm.emapper.annotations                                                  (eggNOG-mapper style)

EVERYTHING HERE IS FAKE (random sequences, invented expression). A database built from these files
is flagged as demo data and the app shows a banner. Standard library only.

    python3 make_demo_data.py --samples ../sample_metadata.csv --out demo_input
"""
import argparse
import csv
import gzip
import json
import math
import os
import random

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TABLE = {a + b + c: AAS[i] for i, (a, b, c) in enumerate((x, y, z) for x in BASES for y in BASES for z in BASES)}
SENSE = [c for c, aa in TABLE.items() if aa != "*"]


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(cds):
    return "".join(TABLE[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))


def rnd_seq(rng, n, gc=0.38):
    return "".join(rng.choices("ACGT", weights=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2], k=n))


def wrap(s, n=70):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


# name, symbol (None -> LOC id), eggNOG preferred name, PFAMs, GO terms, expression pattern
NAMED = [
    ("cadherin-1-like", None, "CDH1", "Cadherin,Cadherin_C_term", "GO:0005509,GO:0007156,GO:0016339", "neural"),
    ("cadherin-2", "Cdh2", "CDH2", "Cadherin,Cadherin_C_term", "GO:0005509,GO:0007156", "neural"),
    ("neural cadherin-like", None, "CDH2", "Cadherin,Cadherin_2", "GO:0005509,GO:0007156", "tailbud"),
    ("protocadherin Fat 4-like", None, "FAT4", "Cadherin,EGF_CA,Laminin_G_2", "GO:0005509,GO:0007156", "gastrula"),
    ("protocadherin-15-like", None, "PCDH15", "Cadherin,Cadherin_2", "GO:0005509", "late"),
    ("cadherin EGF LAG seven-pass G-type receptor 1-like", None, "CELSR1", "Cadherin,7tm_2,Laminin_G_2,EGF_CA", "GO:0004930,GO:0007156", "neural"),
    ("cadherin-related family member 2-like", None, "CDHR2", "Cadherin,Cadherin_3", "GO:0005509", "late"),
    ("desmoglein-2-like", None, "DSG2", "Cadherin,Cadherin_prodomain", "GO:0005509,GO:0007155", "maternal"),
    ("T-box transcription factor TBX2-like", None, "TBX2", "T-box", "GO:0003700,GO:0006355", "gastrula"),
    ("brachyury protein homolog", "Bra", "T", "T-box", "GO:0003700,GO:0006355,GO:0001707", "gastrula"),
    ("homeobox protein Hox-A13-like", "Hox13", "HOXA13", "Homeobox,Hox9_act", "GO:0003700,GO:0009952", "tailbud"),
    ("forkhead box protein A2-like", None, "FOXA2", "Forkhead", "GO:0003700,GO:0006355", "gastrula"),
    ("zinc finger protein 91-like", None, "ZNF91", "zf-C2H2,KRAB", "GO:0003676,GO:0046872", "low"),
    ("transcription factor SOX-9-like", None, "SOX9", "HMG_box", "GO:0003700,GO:0006355", "neural"),
    ("paired box protein Pax-3-like", None, "PAX3", "PAX,Homeobox", "GO:0003700,GO:0007275", "neural"),
    ("GATA-binding factor 4-like", None, "GATA4", "GATA,zf-C2HC", "GO:0003700,GO:0006355", "late"),
    ("fibroblast growth factor 8-like", None, "FGF8", "FGF", "GO:0008083,GO:0008543", "gastrula"),
    ("neurogenic locus notch homolog protein 1-like", None, "NOTCH1", "Notch,EGF,Ank_2", "GO:0007219,GO:0005509", "neural"),
    ("bone morphogenetic protein 2/4-like", None, "BMP2", "TGF_beta,TGFb_propeptide", "GO:0008083,GO:0030509", "gastrula"),
    ("protein Wnt-5a-like", None, "WNT5A", "Wnt", "GO:0005102,GO:0016055", "tailbud"),
    ("frizzled-5-like", None, "FZD5", "Frizzled,Fz", "GO:0004930,GO:0016055", "neural"),
    ("actin, cytoplasmic 1-like", None, "ACTB", "Actin", "GO:0005524,GO:0005856", "ubiquitous"),
    ("tubulin alpha-1A chain-like", None, "TUBA1A", "Tubulin,Tubulin_C", "GO:0003924,GO:0005200", "ubiquitous"),
    ("myosin heavy chain, striated muscle-like", None, "MYH1", "Myosin_head,Myosin_tail_1", "GO:0003774,GO:0016459", "late"),
    ("collagen alpha-1(IV) chain-like", None, "COL4A1", "Collagen,C4", "GO:0005201,GO:0005581", "tailbud"),
    ("histone H3.3-like", None, "H3-3A", "Histone,CENP-T_C", "GO:0003677,GO:0000786", "maternal"),
    ("elongation factor 1-alpha-like", None, "EEF1A1", "GTP_EFTU,GTP_EFTU_D2", "GO:0003746,GO:0005525", "ubiquitous"),
    ("ribosomal protein S3-like", None, "RPS3", "Ribosomal_S3_C,KH_2", "GO:0003735,GO:0006412", "ubiquitous"),
]
GENERIC_PFAM = ["zf-C2H2", "RRM_1", "Pkinase", "WD40", "Ig_3", "EGF", "Ank_2", "LRR_8", "PH", "SH3_1", "Helicase_C", "ABC_tran"]
PATTERNS = {  # relative expression at the 15 stages (egg ... hatching larva)
    "maternal":   [9, 9, 8, 7, 5, 3, 1, 0.6, 0.4, 0.3, 0.2, 0.2, 0.2, 0.1, 0.1],
    "gastrula":   [0.1, 0.1, 0.3, 0.8, 1.5, 3, 8, 9, 5, 3, 2, 1.5, 1, 0.8, 0.5],
    "neural":     [0.1, 0.1, 0.1, 0.2, 0.3, 0.5, 1, 2, 5, 8, 9, 7, 5, 4, 3],
    "tailbud":    [0.1, 0.1, 0.1, 0.1, 0.2, 0.3, 0.5, 0.8, 1.5, 3, 6, 9, 9, 7, 5],
    "late":       [0.05, 0.05, 0.05, 0.1, 0.1, 0.2, 0.3, 0.5, 1, 2, 4, 6, 8, 9, 9],
    "ubiquitous": [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
    "low":        [0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1, 1, 1],
}
MT_GENES = ["COX1", "COX2", "COX3", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "ATP6", "ATP8"]


def make_gene(rng, coding, gene_idx):
    """Return a dict with the plus-strand region sequence and relative (0-based) exon / CDS coordinates."""
    strand = rng.choice("+-")
    n_exons = rng.randint(1, 6)
    if coding:
        aa_len = rng.randint(80, 700)
        cds = "ATG" + "".join(rng.choices(SENSE, k=aa_len - 1)) + rng.choice(["TAA", "TAG", "TGA"])
        utr5, utr3 = rnd_seq(rng, rng.randint(40, 300)), rnd_seq(rng, rng.randint(120, 600))
        tseq = utr5 + cds + utr3
        cds_span = (len(utr5), len(utr5) + len(cds))
    else:
        tseq = rnd_seq(rng, rng.randint(500, 2500))
        cds, cds_span = None, None
    cuts = sorted(rng.sample(range(30, len(tseq) - 30), n_exons - 1)) if n_exons > 1 else []
    bounds = [0] + cuts + [len(tseq)]
    exons_t = [(bounds[i], bounds[i + 1]) for i in range(n_exons)]          # transcript coordinates
    concat, c_off = "", []
    for i, (a, b) in enumerate(exons_t):
        c_off.append(len(concat))
        concat += tseq[a:b]
        if i < n_exons - 1:
            concat += "GT" + rnd_seq(rng, rng.randint(80, 1200)) + "AG"
    region_len = len(concat)

    def to_plus(c0, c1):  # concat coords [c0,c1) -> plus-strand relative coords
        return (c0, c1) if strand == "+" else (region_len - c1, region_len - c0)

    def exon_list(t_end):
        out = []
        for (a, b), c0 in zip(exons_t, c_off):
            if a >= t_end:
                continue
            out.append((a, min(b, t_end), c0))
        return out

    isoforms = []
    iso_ends = [len(tseq)]
    if coding and n_exons >= 2 and rng.random() < 0.5:
        cut = min((len(tseq) - cds_span[1]) // 2, exons_t[-1][1] - exons_t[-1][0] - 5)
        if cut > 10:
            iso_ends.append(len(tseq) - cut)
    for t_end in iso_ends:
        ex = exon_list(t_end)
        exons_rel = [to_plus(c0, c0 + (b - a)) for a, b, c0 in ex]
        cds_rel = []
        if cds_span:
            for a, b, c0 in ex:
                lo, hi = max(a, cds_span[0]), min(b, cds_span[1])
                if lo < hi:
                    cds_rel.append(to_plus(c0 + (lo - a), c0 + (hi - a)))
        isoforms.append({"tseq": tseq[:t_end], "exons": exons_rel, "cds": cds_rel,
                         "cds_seq": cds, "protein": translate(cds)[:-1] if cds else None})
    region = concat if strand == "+" else revcomp(concat)
    return {"strand": strand, "region": region, "isoforms": isoforms}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="sample_metadata.csv")
    ap.add_argument("--out", default="demo_input")
    ap.add_argument("--genes", type=int, default=240)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(a.out, exist_ok=True)

    samples = list(csv.DictReader(open(a.samples)))
    samples.sort(key=lambda r: (int(r["dev_order"]), int(r["replicate"])))

    chroms = {f"NC_9000{i:02d}.1": {"parts": [], "len": 0} for i in range(1, 4)}
    chrom_names = list(chroms)
    gtf, rna, cdsf, prot, eggnog, genes_out = [], [], [], [], [], []
    truth = {"dropped_rna": {}}
    n_named = len(NAMED)
    tx_counter, geneid_counter = 900000, 100000

    def place(chrom, region):
        c = chroms[chrom]
        gap = rnd_seq(rng, rng.randint(300, 2500))
        c["parts"].append(gap)
        c["len"] += len(gap)
        start = c["len"]  # 0-based start of the region
        c["parts"].append(region)
        c["len"] += len(region)
        return start

    for gi in range(a.genes):
        named = NAMED[gi] if gi < n_named else None
        coding = named is not None or rng.random() < 0.88
        g = make_gene(rng, coding, gi)
        chrom = chrom_names[gi % 3]
        start0 = place(chrom, g["region"])
        geneid_counter += 1
        loc = f"LOC{geneid_counter}"
        if named:
            desc, symbol, egg_name, pfams, gos, pattern = named
            gene_id = symbol or loc
        else:
            desc = f"uncharacterized protein {loc}" if coding else f"uncharacterized {loc}"
            gene_id, egg_name, gos = loc, "-", "-"
            pfams = rng.choice(GENERIC_PFAM) if rng.random() < 0.5 else "-"
            pattern = rng.choice(list(PATTERNS))
        biotype = "protein_coding" if coding else "lncRNA"
        gstart = start0 + 1
        gend = start0 + len(g["region"])
        gattr = (f'gene_id "{gene_id}"; transcript_id ""; db_xref "GeneID:{geneid_counter}"; description "{desc}"; '
                 f'gbkey "Gene"; gene "{gene_id}"; gene_biotype "{biotype}";')
        gtf.append(f"{chrom}\tGnomon\tgene\t{gstart}\t{gend}\t.\t{g['strand']}\t.\t{gattr}")
        for iso_i, iso in enumerate(g["isoforms"]):
            tx_counter += 1
            tid = f"{'XM' if coding else 'XR'}_{tx_counter}.1"
            pid = f"XP_{tx_counter}.1" if coding else None
            common = f'gene_id "{gene_id}"; transcript_id "{tid}"; db_xref "GeneID:{geneid_counter}";'
            prod = desc
            tx_s = min(s for s, e in iso["exons"]) + start0 + 1
            tx_e = max(e for s, e in iso["exons"]) + start0
            gtf.append(f"{chrom}\tGnomon\ttranscript\t{tx_s}\t{tx_e}\t.\t{g['strand']}\t.\t{common} gbkey \"mRNA\"; gene "
                       f"\"{gene_id}\"; product \"{prod}\"; transcript_biotype \"{'mRNA' if coding else 'lnc_RNA'}\";")
            order = sorted(range(len(iso["exons"])), key=lambda k: iso["exons"][k][0], reverse=(g["strand"] == "-"))
            # transcription order for exon_number; minus-strand exons are listed in descending order, like NCBI
            for n, k in enumerate(order, 1):
                s, e = iso["exons"][k]
                gtf.append(f"{chrom}\tGnomon\texon\t{s + start0 + 1}\t{e + start0}\t.\t{g['strand']}\t.\t{common} gene "
                           f"\"{gene_id}\"; product \"{prod}\"; exon_number \"{n}\";")
            if coding:
                cds_sorted = sorted(iso["cds"], key=lambda x: x[0], reverse=(g["strand"] == "-"))
                for n, (s, e) in enumerate(cds_sorted, 1):
                    gtf.append(f"{chrom}\tGnomon\tCDS\t{s + start0 + 1}\t{e + start0}\t.\t{g['strand']}\t0\t{common} gene "
                               f"\"{gene_id}\"; product \"{prod}\"; protein_id \"{pid}\"; exon_number \"{n}\";")
            # sequence files; a few transcripts are left out of rna.fna to test the rebuild-from-exons fallback
            if coding and gi in (3, 4) and iso_i == 0:
                truth["dropped_rna"][tid] = iso["tseq"]
            else:
                kind = "PREDICTED: Phallusia mammillata " + desc + f" ({gene_id}), " + ("mRNA" if coding else "ncRNA")
                rna.append(f">{tid} {kind}\n{wrap(iso['tseq'])}")
            if coding:
                loc_str = "join(" + ",".join(f"{s + start0 + 1}..{e + start0}" for s, e in sorted(iso["cds"])) + ")"
                if g["strand"] == "-":
                    loc_str = "complement(" + loc_str + ")"
                cdsf.append(f">lcl|{chrom}_cds_{pid}_{iso_i + 1} [gene={gene_id}] [locus_tag={loc}] [db_xref=GeneID:{geneid_counter}] "
                            f"[protein={desc}] [protein_id={pid}] [location={loc_str}] [gbkey=CDS]\n{wrap(iso['cds_seq'])}")
                prot.append(f">{pid} {desc} [Phallusia mammillata]\n{wrap(iso['protein'])}")
                if (named or rng.random() < 0.6):
                    og = f"COG{rng.randint(1000, 9999)}@1|root,KOG{rng.randint(1000, 9999)}@2759|Eukaryota"
                    eggnog.append("\t".join([
                        pid, f"7719.ENSCING{rng.randint(10**10, 10**11 - 1)}", f"{rng.uniform(1e-120, 1e-8):.2e}",
                        f"{rng.uniform(60, 900):.1f}", og, "Eukaryota", rng.choice("TKSLJOCAUZ"),
                        desc if named else "-", egg_name, gos, "-",
                        f"ko:K{rng.randint(10000, 20000)}" if named else "-", "-", "-", "-", "-", "-", "-", "-", "-", pfams]))
        genes_out.append((gene_id, pattern, coding))

    # mitochondrial "chromosome": exon-only genes like the ones mt_genbank_to_gtf.py writes
    mt_id = "OZ900099.1"
    mt_seq = rnd_seq(rng, 2600, gc=0.45)
    for k, name in enumerate(MT_GENES):
        s, e = 100 + k * 190, 100 + k * 190 + 150
        strand = "+" if k % 3 else "-"
        a_ = f'gene_id "MT-{name}"; gene_name "{name}"; gene_biotype "protein_coding";'
        b_ = f'gene_id "MT-{name}"; transcript_id "MT-{name}.1"; gene_name "{name}"; gene_biotype "protein_coding";'
        gtf.append(f"{mt_id}\tGenBank_MT\tgene\t{s}\t{e}\t.\t{strand}\t.\t{a_}")
        gtf.append(f"{mt_id}\tGenBank_MT\ttranscript\t{s}\t{e}\t.\t{strand}\t.\t{b_}")
        gtf.append(f"{mt_id}\tGenBank_MT\texon\t{s}\t{e}\t.\t{strand}\t.\t{b_}")
        genes_out.append((f"MT-{name}", "ubiquitous", True))
        seq = mt_seq[s - 1:e]
        truth.setdefault("mt", {})[f"MT-{name}.1"] = seq if strand == "+" else revcomp(seq)

    def w(path, text):
        with gzip.open(os.path.join(a.out, path), "wt") as fh:
            fh.write(text)

    w("genome.fna.gz", "".join(f">{n} Phallusia mammillata DEMO chromosome, simulated\n{wrap(''.join(c['parts']))}\n" for n, c in chroms.items())
      + f">{mt_id} Phallusia mammillata DEMO mitochondrion, simulated\n{wrap(mt_seq)}\n")
    w("annotation.gtf.gz", "#gtf-version 2.2\n#DEMO: simulated annotation\n" + "\n".join(gtf) + "\n")
    w("rna.fna.gz", "\n".join(rna) + "\n")
    w("cds.fna.gz", "\n".join(cdsf) + "\n")
    w("proteins.faa.gz", "\n".join(prot) + "\n")

    with open(os.path.join(a.out, "phmamm.emapper.annotations"), "w") as fh:
        fh.write("## emapper-2.1.15 (DEMO, simulated)\n## Tue Oct 27 12:00:00 2026\n")
        fh.write("#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\t"
                 "KEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n")
        fh.write("\n".join(eggnog) + "\n")
        fh.write(f"## {len(prot)} queries scanned\n## Total time (seconds): 1.0\n## Rate: 1 q/s\n")

    # merged-table style gene x sample tables: gene_id, gene_name, then one column per sample
    names = [r["sample"] for r in samples]
    stage_idx = [int(r["dev_order"]) - 1 for r in samples]
    with open(os.path.join(a.out, "salmon.merged.gene_tpm.tsv"), "w") as ft, open(os.path.join(a.out, "salmon.merged.gene_counts.tsv"), "w") as fc:
        hdr = "gene_id\tgene_name\t" + "\t".join(names) + "\n"
        ft.write(hdr)
        fc.write(hdr)
        for gene_id, pattern, coding in genes_out:
            scale = math.exp(rng.gauss(2.2, 1.0))
            tpm_row, cnt_row = [], []
            for si in stage_idx:
                v = scale * PATTERNS[pattern][si] * math.exp(rng.gauss(0, 0.12)) if rng.random() > 0.02 else 0.0
                tpm_row.append(f"{v:.4f}")
                cnt_row.append(f"{max(0.0, v * rng.uniform(15, 25)):.3f}")
            ft.write(f"{gene_id}\t{gene_id}\t" + "\t".join(tpm_row) + "\n")
            fc.write(f"{gene_id}\t{gene_id}\t" + "\t".join(cnt_row) + "\n")

    json.dump(truth, open(os.path.join(a.out, "demo_truth.json"), "w"))
    print(f"DEMO input written to {a.out}: {len(genes_out)} genes, {len(prot)} proteins, {len(eggnog)} eggNOG rows, "
          f"{len(samples)} samples; {len(truth['dropped_rna'])} transcripts left out of rna.fna on purpose")


if __name__ == "__main__":
    main()
