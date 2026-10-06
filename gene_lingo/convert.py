"""Convert gene identifiers between nomenclatures using MyGene.info."""
from . import http
from .species import resolve

MYGENE = "https://mygene.info/v3"
BATCH = 500

# our id-type name -> MyGene field (used for output and as query scope)
ID_TYPES = {
    "symbol": "symbol",
    "name": "name",
    "entrez": "entrezgene",
    "ensembl": "ensembl.gene",
    "hgnc": "HGNC",
    "uniprot": "uniprot.Swiss-Prot",
    "refseq_rna": "refseq.rna",
    "refseq_protein": "refseq.protein",
}
# scopes tried when the input type is not given ("auto"); alias covers old symbols
_AUTO_SCOPES = "symbol,alias,entrezgene,ensembl.gene,HGNC,uniprot,refseq,name"


def _dig(hit, path):
    """Fetch a dotted path from a MyGene hit, flattening lists of dicts; return list of str."""
    cur = [hit]
    for part in path.split("."):
        nxt = []
        for c in cur:
            for item in (c if isinstance(c, list) else [c]):
                if isinstance(item, dict) and part in item:
                    v = item[part]
                    nxt.extend(v if isinstance(v, list) else [v])
        cur = nxt
    return sorted({str(v) for v in cur}, key=str)


def convert(ids, to="symbol", from_type="auto", species="human", fields=None):
    """Map `ids` to identifier type `to`.

    Returns {input_id: [converted ids...]} (empty list if not found).
    `to` may be a list of types, in which case values are {type: [ids]} dicts.
    `from_type` is one of ID_TYPES or "auto". Old/alias symbols resolve via MyGene aliases.
    """
    ids = list(dict.fromkeys(str(i).strip() for i in ids))
    targets = [to] if isinstance(to, str) else list(to)
    for t in targets + ([] if from_type == "auto" else [from_type]):
        if t not in ID_TYPES:
            raise ValueError(f"Unknown id type {t!r}; choose from {sorted(ID_TYPES)}")
    taxid, _ = resolve(species)
    scopes = _AUTO_SCOPES if from_type == "auto" else ID_TYPES[from_type]
    out_fields = ",".join(["symbol"] + [ID_TYPES[t] for t in targets])

    result = {i: {t: [] for t in targets} for i in ids}
    for start in range(0, len(ids), BATCH):
        chunk = ids[start:start + BATCH]
        hits = http.request_json(f"{MYGENE}/query", data={
            "q": ",".join(chunk), "scopes": scopes, "fields": out_fields,
            "species": taxid, "size": 5,
        }) or []
        for hit in hits:
            if hit.get("notfound"):
                continue
            q = hit["query"]
            for t in targets:
                for v in _dig(hit, ID_TYPES[t]):
                    if v not in result[q][t]:
                        result[q][t].append(v)
    if isinstance(to, str):
        return {i: v[to] for i, v in result.items()}
    return result
