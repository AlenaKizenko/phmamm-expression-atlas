#!/usr/bin/env python3
"""Self-test of the database builder: makes simulated inputs, builds a database from them and checks it.

    python3 test_build.py            (standard library only, takes a few seconds)

It checks the parts that are easy to get wrong: strand handling, CDS <-> protein consistency, sequences rebuilt from
exons, identifiers and the search table. It tests the builder against the simulated file formats; the first build
from the real RefSeq / Salmon / eggNOG files should still be checked by eye (the builder prints warnings).
"""
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
B, A = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TABLE = {x + y + z: A[i] for i, (x, y, z) in enumerate((x, y, z) for x in B for y in B for z in B)}
translate = lambda c: "".join(TABLE[c[i:i + 3]] for i in range(0, len(c) - 2, 3))
unz = lambda b: zlib.decompress(b).decode() if b else ""
failed = []


def check(name, ok, extra=""):
    print(("PASS  " if ok else "FAIL  ") + name, extra)
    if not ok:
        failed.append(name)


with tempfile.TemporaryDirectory() as d:
    samples = os.path.join(HERE, "..", "sample_metadata.csv")
    subprocess.run([sys.executable, os.path.join(HERE, "make_demo_data.py"), "--samples", samples, "--out", d], check=True, stdout=subprocess.DEVNULL)
    db_path = os.path.join(d, "test.sqlite")
    subprocess.run([sys.executable, os.path.join(HERE, "build_database.py"), "--gtf", f"{d}/annotation.gtf.gz", "--genome", f"{d}/genome.fna.gz",
                    "--rna", f"{d}/rna.fna.gz", "--cds", f"{d}/cds.fna.gz", "--protein", f"{d}/proteins.faa.gz",
                    "--tpm", f"{d}/salmon.merged.gene_tpm.tsv", "--counts", f"{d}/salmon.merged.gene_counts.tsv", "--samples", samples,
                    "--eggnog", f"{d}/phmamm.emapper.annotations", "--out", db_path, "--demo"], check=True, stdout=subprocess.DEVNULL)
    truth = json.load(open(f"{d}/demo_truth.json"))
    db = sqlite3.connect(db_path)
    q = lambda s, *a: db.execute(s, a).fetchall()

    n_genes, n_tx = q("select count(*) from genes")[0][0], q("select count(*) from transcripts")[0][0]
    check("every gene has a gene sequence, an expression profile (30 samples) and a primary transcript",
          q("select count(*) from seq_gene")[0][0] == n_genes and
          q("select count(*) from (select gene_id from expression group by gene_id having count(*)=30)")[0][0] == n_genes and
          q("select count(*) from genes where primary_transcript is null or primary_transcript=''")[0][0] == 0)
    check("one primary transcript per gene", q("select count(*) from (select gene_id from transcripts where is_primary=1 group by gene_id having count(*)=1)")[0][0] == n_genes)
    rows = q("select t.transcript_id, s.mrna_z, s.cds_z, s.protein from transcripts t join seq_tx s using(transcript_id) where t.protein_id!=''")
    check("CDS translates to the protein (all coding transcripts)", all(translate(unz(c))[:-1] == p for _, _, c, p in rows), f"{len(rows)} transcripts")
    check("CDS lies inside the mRNA (all coding transcripts)", all(unz(c) in unz(m) for _, m, c, _ in rows))
    check("transcripts left out of the RNA FASTA are rebuilt from exons identically",
          all(unz(q("select mrna_z from seq_tx where transcript_id=?", t)[0][0]) == s for t, s in truth["dropped_rna"].items()), f"{len(truth['dropped_rna'])} tested")
    check("mitochondrial (exon-only) transcripts match the simulated sequence, both strands",
          all(unz(q("select mrna_z from seq_tx where transcript_id=?", t)[0][0]) == s for t, s in truth["mt"].items()), f"{len(truth['mt'])} tested")
    check("gene sequence length equals the gene span", all(a == b for _, a, b in q("select g.gene_id, g.end-g.start+1, s.length from genes g join seq_gene s using(gene_id)")))
    plus = q("select gene_id, primary_transcript from genes where strand='+' and biotype='protein_coding' and gene_id not like 'MT-%' limit 1")[0]
    minus = q("select gene_id, primary_transcript from genes where strand='-' and biotype='protein_coding' and gene_id not like 'MT-%' limit 1")[0]
    for label, (gid, tid) in (("+", plus), ("-", minus)):
        gs = unz(q("select seq_zlib from seq_gene where gene_id=?", gid)[0][0])
        m = unz(q("select mrna_z from seq_tx where transcript_id=?", tid)[0][0])
        check(f"{label}-strand gene sequence runs 5'->3' (starts with the mRNA start)", gs.startswith(m[:50]), gid)
    check("exact protein / transcript ID lookups work", q("select gene_id from identifiers where identifier_l='xp_900001.1'") == [("LOC100001",)] and
          q("select gene_id from identifiers where identifier_l='xm_900001'") == [("LOC100001",)])
    cad = {r[0] for r in q("select gene_id from search where symbol_l||' '||desc_l||' '||other_l like '%cadherin%'")}
    check("'cadherin' finds the 8 simulated cadherin-family genes (incl. protocadherins)", len(cad) == 8, str(sorted(cad)))
    check("eggNOG rows are collapsed to one per gene", q("select count(*) from eggnog")[0][0] == q("select count(distinct gene_id) from eggnog")[0][0] > 0)
    check("database is marked as demo", q("select value from meta where key='is_demo'") == [("1",)])

print("\nRESULT:", "all passed" if not failed else f"{len(failed)} FAILED: {failed}")
sys.exit(1 if failed else 0)
