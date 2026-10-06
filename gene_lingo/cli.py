import argparse
import csv
import sys

from . import convert, orthologs, ID_TYPES


def _read_genes(args):
    genes = list(args.genes)
    if args.file:
        with open(args.file) as fh:
            genes += [l.strip() for l in fh if l.strip()]
    if not genes:
        genes = [l.strip() for l in sys.stdin if l.strip()]
    return genes


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

    a = p.parse_args(argv)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
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


if __name__ == "__main__":
    main()
