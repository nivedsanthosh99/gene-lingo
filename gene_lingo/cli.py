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
    m = bacteria.Mapper(a.tables)
    dst = m.key_for(a.dst)
    w.writerow(["input", "locus_tag", m.genomes[dst]["label"], "pident"])
    for r in m.map(_read_genes(a), a.src, [dst]):
        w.writerow([r["input"], r["locus_tag"], r[dst], "" if r[f"{dst}_pident"] is None else r[f"{dst}_pident"]])


def _build_args(p):
    p.add_argument("-o", "--outdir", default="bacteria_tables")
    p.add_argument("--cache", help="download cache (default ~/.cache/gene-lingo)")
    p.add_argument("--threads", type=int, default=4)


def _progress(msg):
    print(msg, file=sys.stderr)


def _summarise(prot, hub_map, registry, outdir, hub):
    from . import bacteria
    print(f"hub: {registry[hub]['label']}; wrote {outdir}/ (hub.tsv, core.tsv, pairs/, genes.tsv.gz, genomes.json)", file=sys.stderr)
    for k, m in hub_map.items():
        print(f"  {registry[k]['label']}: {len(m)} of {len(prot[hub].cds)} hub genes have an ortholog", file=sys.stderr)
    print(f"  core across {len({g['species'] for g in registry.values()})} species: "
          f"{len(bacteria.species_core(hub_map, registry))} genes", file=sys.stderr)


def _genomes(a, w):
    from . import bacteria, genomes
    cache = a.cache if getattr(a, "cache", None) else genomes.default_cache()
    if a.gcmd == "find":
        w.writerow(["accession", "organism", "strain", "category", "level", "taxid"])
        for r in genomes.find_assemblies(a.query or None, a.taxid, cache, a.n):
            w.writerow([r["acc"], r["organism"], r["strain"], r["category"], r["level"], r["taxid"]])
        return
    idx = genomes.assembly_index(cache)
    registry = dict(genomes.parse_genome_arg(x, cache, idx) for x in a.specs)
    hub = None
    if a.hub:
        hub, spec = genomes.parse_genome_arg(a.hub, cache, idx)
        registry.setdefault(hub, spec)
    hub = hub or next(iter(registry))
    prot, hub_map, _ = bacteria.build_set(registry, a.outdir, os.path.join(cache, "genomes"), hub, a.hub_only, a.threads, _progress)
    _summarise(prot, hub_map, registry, a.outdir, hub)


def _osdr(a, w):
    from . import bacteria, genomes, osdr
    cache = a.cache if getattr(a, "cache", None) else genomes.default_cache()
    if a.ocmd == "survey":
        rows = osdr.survey(cache, os.path.join(cache, "genomes"), progress=_progress)
        osdr.write_survey(rows, a.out)
        print(f"wrote {len(rows)} studies to {a.out}", file=sys.stderr)
        return
    if a.ocmd == "info":
        registry, picks = osdr.genomes_for(a.studies, cache, os.path.join(cache, "genomes"), _progress)
        w.writerow(["osd", "organism", "strain", "processed_table", "genome", "label", "confidence", "evidence"])
        for n, p in picks.items():
            i = p["info"]
            w.writerow([f"OSD-{n}", "; ".join(i["organisms"]), "; ".join(i["strains"]),
                        i["tables"][0][0] if i["tables"] else "(raw data only)", p["acc"], p["label"], p["confidence"], p["evidence"]])
        return
    registry, picks = osdr.genomes_for(a.studies, cache, os.path.join(cache, "genomes"), _progress)
    for n, p in picks.items():
        print(f"OSD-{n}: {p['label']} ({p['acc']}) [{p['confidence']}] {p['evidence']}", file=sys.stderr)
    hub = next((k for k, g in registry.items() if g["species"] == "Bacillus subtilis"), next(iter(registry)))
    prot, hub_map, _ = bacteria.build_set(registry, a.outdir, os.path.join(cache, "genomes"), hub, a.hub_only, a.threads, _progress)
    _summarise(prot, hub_map, registry, a.outdir, hub)


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
    bm.add_argument("--from", dest="src", required=True, help="genome key, label or OSD study, e.g. OSD-554")
    bm.add_argument("--to", dest="dst", required=True)
    bm.add_argument("-d", "--tables", default="bacteria_tables")

    gn = sub.add_parser("genomes", help="orthologs/locus-tag mapping for any bacterial genomes")
    gsub = gn.add_subparsers(dest="gcmd", required=True)
    gf = gsub.add_parser("find", help="search NCBI RefSeq for a genome by organism/strain name or taxid")
    gf.add_argument("query", nargs="?", default="")
    gf.add_argument("--taxid")
    gf.add_argument("-n", type=int, default=5)
    gb = gsub.add_parser("build", help="map genes between genomes (accession, 'Genus species strain', taxid:N, or gff,faa files)")
    gb.add_argument("specs", nargs="+")
    gb.add_argument("--hub", help="genome spec to use as the hub (default: the first)")
    gb.add_argument("--hub-only", action="store_true", help="only compare each genome with the hub (faster for many genomes)")
    _build_args(gb)

    os_ = sub.add_parser("osdr", help="resolve and map the genomes of any OSDR study")
    osub = os_.add_subparsers(dest="ocmd", required=True)
    oi = osub.add_parser("info", help="organism, strain, processed tables and chosen genome of a study")
    oi.add_argument("studies", nargs="+")
    ov = osub.add_parser("survey", help="genome choice for every bacterial transcription study in OSDR")
    ov.add_argument("-o", "--out", default="osdr_bacterial_survey.tsv")
    ov.add_argument("--cache")
    ob = osub.add_parser("build", help="map genes between the genomes of several studies")
    ob.add_argument("studies", nargs="+")
    ob.add_argument("--hub-only", action="store_true")
    _build_args(ob)


    a = p.parse_args(argv)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    if a.cmd == "bacteria":
        return _bacteria(a, w)
    if a.cmd == "genomes":
        return _genomes(a, w)
    if a.cmd == "osdr":
        return _osdr(a, w)
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
    except (urllib.error.URLError, ValueError, KeyError, RuntimeError) as e:
        sys.exit(f"gene-lingo: {e}")


if __name__ == "__main__":
    run()
