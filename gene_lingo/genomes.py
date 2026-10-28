"""Find and describe bacterial/archaeal genomes for any organism (not tied to OSDR).

A genome spec is what `bacteria.build(genomes=...)` takes: {key: {label, species, acc, osds, [gff, faa]}}.
`parse_genome_arg` turns user input into one: an NCBI assembly accession, an organism/strain name,
`taxid:<n>`, or a pair of local files `genes.gff.gz,proteins.faa.gz[,label]`.
"""
import csv
import gzip
import os
import re
import urllib.request

REFSEQ_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
ACC_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")
CATEGORY_RANK = {"reference genome": 3, "representative genome": 2}
COLS = ["acc", "category", "taxid", "species_taxid", "organism", "strain", "level", "asm_name", "date", "group"]


def default_cache():
    return os.environ.get("GENE_LINGO_CACHE") or os.path.join(os.path.expanduser("~"), ".cache", "gene-lingo")


def assembly_index(cache=None, url=REFSEQ_SUMMARY, force=False):
    """Slim local copy of NCBI's RefSeq assembly summary: current bacterial/archaeal complete genomes plus
    all reference/representative genomes.  Downloaded once (~250 MB, streamed) and cached."""
    cache = cache or default_cache()
    path = os.path.join(cache, "assembly_index.tsv.gz")
    if force or not os.path.exists(path):
        os.makedirs(cache, exist_ok=True)
        tmp = path + ".part"
        with urllib.request.urlopen(url, timeout=120) as resp, gzip.open(tmp, "wt") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(COLS)
            for raw in resp:
                line = raw.decode("utf-8", "replace")
                if line.startswith("#") or line.count("\t") < 20:
                    continue
                f = line.rstrip("\n").split("\t")
                if f[10] != "latest" or f[24] not in ("bacteria", "archaea"):
                    continue
                if f[11] not in ("Complete Genome", "Chromosome") and f[4] not in CATEGORY_RANK:
                    continue
                strain = re.sub(r"^(strain|isolate|serovar|cultivar)=", "", f[8]) or f[9]
                w.writerow([f[0], f[4], f[5], f[6], f[7], strain, f[11], f[15], f[14], f[24]])
        os.replace(tmp, path)
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def species_name(organism):
    """'Escherichia coli str. K-12 substr. MG1655' -> 'Escherichia coli'."""
    return " ".join(organism.split()[:2])


def _rank(row, strain_hit):
    return (strain_hit, CATEGORY_RANK.get(row["category"], 0), row["level"] == "Complete Genome", row["date"])


def _strain_match(strain, words):
    """3: a query word is the whole strain name; 2: it is one of the strain's words; 0: no match."""
    st = strain.lower()
    if not st or st == "na":
        return 0
    if any(w == st for w in words):
        return 3
    return 2 if any(w in re.split(r"[\s/_]", st) for w in words) else 0


def find_assemblies(query=None, taxid=None, cache=None, limit=10, index=None):
    """Candidate assemblies, best first.

    `taxid` matches the strain-level or species-level taxon; `query` is free text ("Aliivibrio fischeri ES114",
    "UAMS-1"): every word must occur in the organism name / strain, and a strain-word match ranks first.
    """
    rows = index if index is not None else assembly_index(cache)
    hits = []
    if taxid:
        t = str(taxid)
        for r in rows:
            if r["taxid"] == t or r["species_taxid"] == t:
                # a named strain whose own taxon is the study's beats a species-level match
                named = (r["taxid"] == t and r["species_taxid"] != t and r["strain"] not in ("", "na")
                         and r["strain"].lower() in r["organism"].lower())
                hits.append((_rank(r, 2 if named else 1), r))
    if query:
        words = [w for w in re.split(r"[\s_]+", query.lower()) if w and w not in ("strain", "str", "str.", "substr.", "subsp.", "sp.")]
        for r in rows:
            hay = f"{r['organism']} {r['strain']}".lower()
            if all(w in hay for w in words):
                hit = _strain_match(r["strain"], words[2:]) if len(words) > 2 else 1
                hits.append((_rank(r, hit), r))
    seen, out = set(), []
    for _, r in sorted(hits, key=lambda x: x[0], reverse=True):
        if r["acc"] not in seen:
            seen.add(r["acc"])
            out.append(r)
    return out[:limit]


def lookup(acc, cache=None, index=None):
    """Index row for an assembly accession (matching the digits, so GCA_/GCF_ and versions interchange)."""
    rows = index if index is not None else assembly_index(cache)
    digits = re.sub(r"\D", "", acc.split(".")[0])
    best = [r for r in rows if r["acc"][4:13] == digits]
    return sorted(best, key=lambda r: (r["acc"] == acc, r["acc"]), reverse=True)[0] if best else None


def label_for(row):
    """'E. coli K-12 MG1655'-style short label from an index row."""
    sp = species_name(row["organism"]).split()
    base = f"{sp[0][0]}. {sp[1]}" if len(sp) > 1 else row["organism"]
    return f"{base} {row['strain']}".strip() if row["strain"] else base


def spec_from_row(row, osds=()):
    return dict(label=label_for(row), species=species_name(row["organism"]), acc=row["acc"], osds=list(osds))


def parse_genome_arg(arg, cache=None, index=None):
    """-> (key, spec).  See the module docstring for the accepted forms."""
    arg = arg.strip()
    if "," in arg and ".gff" in arg.lower():
        parts = [p.strip() for p in arg.split(",")]
        label = parts[2] if len(parts) > 2 else os.path.basename(parts[0]).split(".")[0]
        key = re.sub(r"\W+", "_", label).strip("_")
        return key, dict(label=label, species=label, acc="", osds=[], gff=parts[0], faa=parts[1])
    if ACC_RE.match(arg):
        row = None
        try:
            row = lookup(arg, cache, index)
        except OSError:
            pass
        spec = spec_from_row(row) if row else dict(label=arg, species=arg, osds=[])
        spec["acc"] = arg
        return arg, spec
    if arg.lower().startswith("taxid:"):
        hits = find_assemblies(taxid=arg.split(":", 1)[1], cache=cache, index=index)
    else:
        hits = find_assemblies(query=arg, cache=cache, index=index)
    if not hits:
        raise ValueError(f"no RefSeq bacterial/archaeal genome found for {arg!r}")
    return hits[0]["acc"], spec_from_row(hits[0])
