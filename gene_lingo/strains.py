"""Between-strain comparison: genomes of one species, from the pair tables built by `bacteria build`.

Strains of one species share almost all genes with near-identical proteins, so the question is different from the
between-species one: which genes match, and which genes exist in only one strain.
"""
import statistics

# identity bands for the histogram: (lower bound, upper bound, label)
BANDS = [(0, 90, "below 90%"), (90, 95, "90-95%"), (95, 99, "95-99%"), (99, 99.9, "99-99.9%"), (99.9, 100.01, "99.9-100%")]


def species_groups(mapper):
    """{species: [genome keys]} for species that have at least two genomes (strains) in the tables."""
    groups = {}
    for key, g in mapper.genomes.items():
        groups.setdefault(g["species"], []).append(key)
    return {sp: ks for sp, ks in groups.items() if len(ks) > 1}


def compare(mapper, a, b):
    """Compare strain genomes a and b.

    -> dict with `shared` [(tag_a, tag_b, pident)], `only_a` / `only_b` (locus tags with no counterpart found),
    gene counts, median identity and the share of shared genes at >= 99 % identity.
    """
    pair = mapper.pair(a, b)
    shared = [(ta, tb, pid) for ta, (tb, pid) in pair.items()]
    matched_b = {tb for _, tb, _ in shared}
    ids = [pid for _, _, pid in shared]
    return dict(
        a=a, b=b, shared=shared, n_a=len(mapper.genes[a].cds), n_b=len(mapper.genes[b].cds),
        only_a=[t for t in mapper.genes[a].cds if t not in pair],
        only_b=[t for t in mapper.genes[b].cds if t not in matched_b],
        median_identity=statistics.median(ids) if ids else None,
        share_ge_99=(sum(x >= 99 for x in ids) / len(ids)) if ids else None)


def identity_histogram(shared):
    """[(band label, number of shared genes)] in BANDS order."""
    return [(label, sum(lo <= pid < hi for _, _, pid in shared)) for lo, hi, label in BANDS]


def all_pairs(mapper):
    """One summary row per pair of strains of the same species."""
    rows = []
    for sp, keys in species_groups(mapper).items():
        for i, a in enumerate(keys):
            for b in keys[i + 1:]:
                c = compare(mapper, a, b)
                rows.append(dict(species=sp, strain_a=mapper.genomes[a]["label"], strain_b=mapper.genomes[b]["label"],
                                 genes_a=c["n_a"], genes_b=c["n_b"], shared=len(c["shared"]), only_a=len(c["only_a"]),
                                 only_b=len(c["only_b"]), median_identity=c["median_identity"], share_ge_99=c["share_ge_99"]))
    return rows
