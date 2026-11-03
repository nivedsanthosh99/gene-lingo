"""Make data/example_alignment.json: a real DIAMOND comparison of one protein, for the Explanation page.

Takes the PA14 DnaA protein (PA14_RS00005), finds its best hits among the B. subtilis 168 proteins, then does
the reverse search from the top hit, exactly as the ortholog step does (very-sensitive, e-value <= 1e-5).
Needs the `diamond` binary and (for the similarity marks only) Biopython:  pip install biopython

    python scripts/make_example_alignment.py [cache_dir]
"""
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from gene_lingo import bacteria, genomes  # noqa: E402

QUERY_GENOME, SUBJECT_GENOME, QUERY_TAG = "GCF_000014625.1", "GCF_000009045.1", "PA14_RS00005"
FMT = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore", "qstart", "qend", "sstart", "send",
       "qseq_gapped", "sseq_gapped"]


def load(acc, cache):
    return bacteria.Proteome(acc, *bacteria.fetch_genome(acc, cache, {acc: dict(acc=acc)}))


def search(query_prot, tag, subject_prot, work):
    qf, sf = os.path.join(work, "q.faa"), os.path.join(work, "s.faa")
    with open(qf, "w") as fh:
        fh.write(f">{tag}\n{query_prot.cds[tag]['seq']}\n")
    subject_prot.write_fasta(sf)
    subprocess.run(["diamond", "makedb", "--in", sf, "-d", os.path.join(work, "s"), "--quiet"], check=True)
    out = os.path.join(work, "out.tsv")
    subprocess.run(["diamond", "blastp", "-q", qf, "-d", os.path.join(work, "s"), "-o", out, "--evalue", "1e-5", "-k", "5",
                    "--very-sensitive", "--quiet", "--outfmt", "6", *FMT], check=True)
    rows = [dict(zip(FMT, line.rstrip("\n").split("\t"))) for line in open(out)]
    return sorted(rows, key=lambda r: -float(r["bitscore"]))


def describe(rows, subject_prot):
    out = []
    for r in rows:
        ql, sl, ln = int(r["qlen"]), int(r["slen"]), int(r["length"])
        rec = subject_prot.cds[r["sseqid"]]
        # shown: the share of each protein's residues that lines up; the pipeline's rule (passes_coverage) divides the
        # alignment length, which also counts gap columns, by each protein's length
        out.append(dict(subject=r["sseqid"], symbol=rec["symbol"], product=rec["product"], pident=float(r["pident"]),
                        aligned=ln, query_len=ql, subject_len=sl,
                        query_coverage=round((int(r["qend"]) - int(r["qstart"]) + 1) / ql, 2),
                        subject_coverage=round((int(r["send"]) - int(r["sstart"]) + 1) / sl, 2),
                        evalue=float(r["evalue"]), bitscore=float(r["bitscore"]),
                        passes_coverage=ln / ql >= 0.5 and ln / sl >= 0.5))
    return out


def marks(q, s):
    """'|' identical, ':' similar (BLOSUM62 > 0), ' ' different, gaps blank."""
    try:
        from Bio.Align import substitution_matrices
        m = substitution_matrices.load("BLOSUM62")
    except ImportError:
        m = None
    line = []
    for a, b in zip(q, s):
        if a == "-" or b == "-":
            line.append(" ")
        elif a == b:
            line.append("|")
        elif m is not None and m[a][b] > 0:
            line.append(":")
        else:
            line.append(" ")
    return "".join(line)


def main():
    cache = (sys.argv[1] if len(sys.argv) > 1 else genomes.default_cache()).rstrip("/")
    cache = cache if cache.endswith("genomes") else os.path.join(cache, "genomes")
    qp, sp = load(QUERY_GENOME, cache), load(SUBJECT_GENOME, cache)
    with tempfile.TemporaryDirectory() as w1, tempfile.TemporaryDirectory() as w2:
        fwd = search(qp, QUERY_TAG, sp, w1)
        best = fwd[0]["sseqid"]
        rev = search(sp, best, qp, w2)
    top = fwd[0]
    qs, ss = top["qseq_gapped"], top["sseq_gapped"]
    mk = marks(qs, ss)
    blocks = []
    qpos, spos = int(top["qstart"]), int(top["sstart"])
    for i in range(0, min(len(qs), 180), 60):
        q, m, s = qs[i:i + 60], mk[i:i + 60], ss[i:i + 60]
        blocks.append(dict(q=q, m=m, s=s, qpos=qpos, spos=spos))
        qpos += len(q) - q.count("-")
        spos += len(s) - s.count("-")
    data = dict(
        query=dict(tag=QUERY_TAG, genome=QUERY_GENOME, symbol=qp.cds[QUERY_TAG]["symbol"], product=qp.cds[QUERY_TAG]["product"],
                   length=len(qp.cds[QUERY_TAG]["seq"])),
        subject_genome=SUBJECT_GENOME,
        params=dict(program="DIAMOND blastp --very-sensitive", evalue=1e-5, min_coverage=0.5),
        forward=describe(fwd, sp), reverse=describe(rev, qp), reverse_query=best,
        reciprocal=bool(rev) and rev[0]["sseqid"] == QUERY_TAG,
        alignment=dict(subject=top["sseqid"], pident=float(top["pident"]), length=int(top["length"]), columns_shown=min(len(qs), 180),
                       qstart=int(top["qstart"]), qend=int(top["qend"]), sstart=int(top["sstart"]), send=int(top["send"]),
                       identical=mk.count("|"), similar=mk.count(":"), gaps=qs.count("-") + ss.count("-"), blocks=blocks))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "example_alignment.json")
    with open(out, "w") as fh:
        json.dump(data, fh, indent=1)
    print("wrote", os.path.normpath(out), "reciprocal:", data["reciprocal"], "forward hits:", len(fwd), "reverse hits:", len(rev))


if __name__ == "__main__":
    main()
