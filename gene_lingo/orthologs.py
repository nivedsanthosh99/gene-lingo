"""Ortholog lookup via the Ensembl Compara REST API."""
from . import http
from .convert import convert
from .species import resolve

ENSEMBL = "https://rest.ensembl.org"


def _is_ensembl(gene):
    return gene.upper().startswith("ENS") and any(c.isdigit() for c in gene)


def orthologs(genes, source="human", target="mouse", ortholog_type="all", symbols=True):
    """Find `target`-species orthologs of `genes` (symbols, Entrez ids or Ensembl ids).

    Returns {input_gene: [ {id, symbol, species, type, perc_id, perc_pos} ... ]}.
    ortholog_type: 'all', 'one2one', 'one2many' or 'many2many'.
    Inputs that are not Ensembl IDs are first converted to Ensembl gene IDs via MyGene.
    """
    _, src_ens = resolve(source)
    _, tgt_ens = resolve(target)
    genes = list(dict.fromkeys(str(g).strip() for g in genes))

    to_lookup = [g for g in genes if not _is_ensembl(g)]
    ens_ids = convert(to_lookup, to="ensembl", species=source) if to_lookup else {}

    result = {}
    for g in genes:
        gids = [g] if _is_ensembl(g) else ens_ids.get(g, [])
        found = []
        for gid in gids:
            data = http.request_json(
                f"{ENSEMBL}/homology/id/{src_ens}/{gid}",
                params={"type": "orthologues", "target_species": tgt_ens,
                        "sequence": "none", "content-type": "application/json"})
            for rec in (data or {}).get("data", []):
                for h in rec.get("homologies", []):
                    if ortholog_type != "all" and h["type"] != f"ortholog_{ortholog_type}":
                        continue
                    t = h["target"]
                    found.append({"id": t["id"], "symbol": None, "species": t["species"],
                                  "type": h["type"], "perc_id": t.get("perc_id"),
                                  "perc_pos": t.get("perc_pos")})
        result[g] = found

    if symbols:
        all_ids = sorted({o["id"] for v in result.values() for o in v})
        sym = convert(all_ids, to="symbol", from_type="ensembl", species=target) if all_ids else {}
        for v in result.values():
            for o in v:
                o["symbol"] = (sym.get(o["id"]) or [None])[0]
    return result
