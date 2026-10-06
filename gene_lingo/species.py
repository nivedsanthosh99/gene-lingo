"""Species name resolution: common name / Latin name / taxid -> (taxid, ensembl name)."""

# common name, taxid, ensembl species name
_SPECIES = [
    ("human", 9606, "homo_sapiens"),
    ("mouse", 10090, "mus_musculus"),
    ("rat", 10116, "rattus_norvegicus"),
    ("zebrafish", 7955, "danio_rerio"),
    ("fly", 7227, "drosophila_melanogaster"),
    ("worm", 6239, "caenorhabditis_elegans"),
    ("yeast", 559292, "saccharomyces_cerevisiae"),
    ("chicken", 9031, "gallus_gallus"),
    ("pig", 9823, "sus_scrofa"),
    ("cow", 9913, "bos_taurus"),
    ("dog", 9615, "canis_lupus_familiaris"),
    ("macaque", 9544, "macaca_mulatta"),
    ("chimpanzee", 9598, "pan_troglodytes"),
    ("frog", 8364, "xenopus_tropicalis"),
]
_ALIASES = {"drosophila": "fly", "c. elegans": "worm", "celegans": "worm",
            "saccharomyces cerevisiae": "yeast", "homo sapiens": "human",
            "mus musculus": "mouse", "danio rerio": "zebrafish",
            "rattus norvegicus": "rat", "cattle": "cow", "xenopus": "frog"}


def resolve(species):
    """Return (taxid, ensembl_name). Accepts 'mouse', 'Mus musculus', 'mus_musculus', 10090."""
    key = str(species).strip().lower().replace("_", " ")
    key = _ALIASES.get(key, key)
    for common, taxid, ens in _SPECIES:
        if key in (common, str(taxid), ens.replace("_", " ")):
            return taxid, ens
    raise ValueError(f"Unknown species {species!r}; known: {', '.join(s[0] for s in _SPECIES)}")
