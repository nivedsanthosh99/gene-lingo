import argparse
import csv
import os
import sys
import tempfile
import urllib.error

from . import convert, orthologs, ID_TYPES


def _read_genes(args):
    genes = list(args.genes)
    if args.file:
        with open(args.file) as fh:
            genes += [l.strip() for l in fh if l.strip()]
    if not genes:
        genes = [l.strip() for l in sys.stdin if l.strip()]
    return genes


def _bacteria(a, w):
    from . import bacteria
    if a.bcmd == "build":
        with tempfile.TemporaryDirectory() as work:
            prot, hub, pair = bacteria.build(a.cache, work, threads=a.threads)
        bacteria.write_tables(a.outdir, prot, hub, pair)
        print(f"wrote {a.outdir}/hub.tsv, core.tsv ({len(bacteria.species_core(hub))} genes), pairs/", file=sys.stderr)
        return
    src, dst = bacteria.genome_for_osd(a.src), bacteria.genome_for_osd(a.dst)
    prot = bacteria.Proteome(src, *bacteria.fetch_genome(src, a.cache))
    pair = {} if src == dst else bacteria.load_pair(a.tables, src, dst)
    w.writerow(["input", f"{dst}_locus_tag"])
    for g, hit in bacteria.translate(_read_genes(a), src, dst, pair, prot).items():
        w.writerow([g, hit or ""])


def main(argv=None):
    p = argparse.ArgumentParser(prog="gene-lingo")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("convert", help="convert gene ids between nomenclatures")
    c.add_argument("genes", nargs="*")
    c.add_argument("-f", "--file", help="file with one gene per line (default: stdin)")
    c.add_argument("--from", dest="from_type", default="auto", choices=["auto", *ID_TYPES])
    c.add_argument("--to", default="symbol", choices=list(ID_TYPES))
    c.add_argument("-s", "--species", default="human")

    o = sub.add_parser("orthologs", help="find orthologs in another species")
    o.add_argument("genes", nargs="*")
    o.add_argument("-f", "--file")
    o.add_argument("-s", "--source", default="human")
    o.add_argument("-t", "--target", default="mouse")
    o.add_argument("--type", default="all", choices=["all", "one2one", "one2many", "many2many"])

    bc = sub.add_parser("bacteria", help="map locus tags between the bacterial genomes used by OSDR studies")
    bsub = bc.add_subparsers(dest="bcmd", required=True)
    bb = bsub.add_parser("build", help="download genomes, compute orthologs, write tables")
    bb.add_argument("-o", "--outdir", default="bacteria_tables")
    bb.add_argument("--cache", default="genomes_cache")
    bb.add_argument("--threads", type=int, default=4)
    bm = bsub.add_parser("map", help="translate locus tags from one OSD study to another")
    bm.add_argument("genes", nargs="*")
    bm.add_argument("-f", "--file")
    bm.add_argument("--from", dest="src", required=True, help="OSD number, e.g. OSD-554")
    bm.add_argument("--to", dest="dst", required=True)
    bm.add_argument("-d", "--tables", default="bacteria_tables")
    bm.add_argument("--cache", default="genomes_cache")

    a = p.parse_args(argv)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    if a.cmd == "bacteria":
        return _bacteria(a, w)
    genes = _read_genes(a)
    if a.cmd == "convert":
        w.writerow(["input", a.to])
        for g, vals in convert(genes, a.to, a.from_type, a.species).items():
            w.writerow([g, ";".join(vals)])
    else:
        w.writerow(["input", "ortholog_id", "ortholog_symbol", "species", "type", "perc_id"])
        for g, hits in orthologs(genes, a.source, a.target, a.type).items():
            for h in hits or [{}]:
                w.writerow([g, h.get("id", ""), h.get("symbol") or "", h.get("species", ""),
                            h.get("type", ""), h.get("perc_id", "")])


def run():
    try:
        main()
    except (urllib.error.URLError, ValueError) as e:
        sys.exit(f"gene-lingo: {e}")


if __name__ == "__main__":
    run()
