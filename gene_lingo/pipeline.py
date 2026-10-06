"""One-call translation: (input species, id type) -> (output species, id type)."""
from .convert import convert
from .orthologs import orthologs
from .species import resolve


def translate(genes, in_species="human", out_species="human", in_type="auto", out_type="symbol",
              ortholog_type="all"):
    """Return a list of row dicts; one row per (input gene, output id), blank output if unmapped.

    Same species: plain ID conversion. Different species: find orthologs, then express
    them in `out_type`.
    """
    genes = list(dict.fromkeys(str(g).strip() for g in genes if str(g).strip()))
    in_taxid, in_ens = resolve(in_species)
    out_taxid, out_ens = resolve(out_species)
    base = {"input_species": in_ens, "output_species": out_ens}
    rows = []

    if in_taxid == out_taxid:
        for g, vals in convert(genes, to=out_type, from_type=in_type, species=in_species).items():
            for v in vals or [""]:
                rows.append({"input": g, **base, "output_id": v, "orthology_type": "", "perc_id": ""})
        return rows

    # to find orthologs we need Ensembl ids in the source species
    hits = orthologs(genes, in_species, out_species, ortholog_type, symbols=False) \
        if in_type in ("auto", "symbol", "ensembl", "entrez", "hgnc") else None
    if hits is None:
        raise ValueError(f"Orthology search needs a gene-level input type, not {in_type!r}")
    ens = sorted({h["id"] for v in hits.values() for h in v})
    conv = convert(ens, to=out_type, from_type="ensembl", species=out_species) \
        if ens and out_type != "ensembl" else {e: [e] for e in ens}
    for g, found in hits.items():
        if not found:
            rows.append({"input": g, **base, "output_id": "", "orthology_type": "", "perc_id": ""})
        for h in found:
            for v in conv.get(h["id"]) or [""]:
                rows.append({"input": g, **base, "output_id": v, "orthology_type": h["type"],
                             "perc_id": h["perc_id"] if h["perc_id"] is not None else ""})
    return rows
