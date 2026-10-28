"""Work out, for any OSDR study, which genome its gene identifiers belong to.

OSDR metadata gives the organism (NCBI taxon), usually the strain, and, for studies GeneLab processed, the
assembly used.  Processed count/probe tables give the actual identifiers.  `describe` gathers that,
`choose_genome` picks the genome (checking identifiers against each candidate when a table exists) and
`survey` does it for every bacterial transcriptomics study.
"""
import csv
import functools
import gzip
import hashlib
import io
import json
import os
import re
import time
import urllib.error
import urllib.request

from . import bacteria, genomes

API = "https://osdr.nasa.gov/osdr/data"
BASE = "https://osdr.nasa.gov"
MAX_TABLE_BYTES = 80_000_000
TIE_MARGIN = 0.15  # a better-named candidate within this much of the best identifier match wins
MIN_MATCH = 0.45  # share of a table's ids a genome must recognise (the rest are rRNA/tRNA/pseudogenes, controls)


@functools.lru_cache(maxsize=16)
def _get(url, retries=3):
    req = urllib.request.Request(url, headers={"User-Agent": "gene-lingo/0.1"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except (urllib.error.URLError, OSError):     # includes timeouts and connection resets
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def get_json(url, cache=None, max_age=7 * 86400):
    """JSON from the OSDR API; with `cache` (a directory) study metadata/file lists are kept on disk for a week."""
    path = None
    if cache and re.search(r"/osd/(meta|files)/", url):
        os.makedirs(cache, exist_ok=True)
        path = os.path.join(cache, hashlib.md5(url.encode()).hexdigest() + ".json")
        if os.path.exists(path) and time.time() - os.path.getmtime(path) < max_age:
            with open(path) as fh:
                return json.load(fh)
    data = json.loads(_get(url))
    if path:
        with open(path, "w") as fh:
            json.dump(data, fh)
    return data


def osd_number(osd):
    m = re.fullmatch(r"\s*(?:OSD-?|GLDS-?)?(\d+)\s*", str(osd), re.I)
    if not m:
        raise ValueError(f"not an OSD study id: {osd!r}")
    return int(m.group(1))


def search_studies(term="", batch=100):
    """Every OSD-numbered study matching `term` (GeneLab source), as the search API's `_source` dicts."""
    out, start = [], 0
    while True:
        d = get_json(f"{API}/search?term={urllib.request.quote(term)}&from={start}&size={batch}&type=cgene")
        hits = d["hits"]["hits"]
        out += [h["_source"] for h in hits]
        start += batch
        if not hits or start >= d["hits"]["total"]:
            return out


# ---- metadata -----------------------------------------------------------------

def _walk(o):
    if isinstance(o, dict):
        yield o
        for v in o.values():
            yield from _walk(v)
    elif isinstance(o, list):
        for v in o:
            yield from _walk(v)


def parse_meta(meta, osd):
    """Organism(s), taxids, strain text, accession hints and processing text from an OSDR meta record."""
    study = meta["study"][f"OSD-{osd}"]["studies"][0]
    taxa = {}
    for d in _walk(study):
        m = re.search(r"NCBITAXON[_/](\d+)", str(d.get("termAccession", "")), re.I)
        if m and d.get("annotationValue"):
            taxa[m.group(1)] = d["annotationValue"]
    cats = {c["@id"]: (c.get("characteristicType") or {}).get("annotationValue", "") for c in study.get("characteristicCategories", [])}
    strains, organisms = set(), set()
    for d in _walk(study.get("materials", {})):
        cid = (d.get("category") or {}).get("@id") if isinstance(d.get("category"), dict) else None
        if cid and re.fullmatch(r"organism", cats.get(cid, ""), re.I):
            v = d.get("value")
            v = v.get("annotationValue") if isinstance(v, dict) else v
            if v and str(v).strip() and str(v).strip().lower() not in ("not applicable", "na", "n/a"):
                organisms.add(str(v).strip())
        if cid and re.search(r"strain|isolate|genotype", cats.get(cid, ""), re.I):
            v = d.get("value")
            v = v.get("annotationValue") if isinstance(v, dict) else v
            if v and str(v).strip() and str(v).lower() not in ("not applicable", "na", "n/a"):
                strains.add(str(v).strip())
    text = json.dumps(study)
    accs = sorted(set(re.findall(r"GC[AF]_\d{9}\.\d+", text)))
    digits = sorted(set(re.findall(r"gc[af]_(\d{9})", text, re.I)) - {a[4:13] for a in accs})
    return dict(title=study.get("title", ""), taxa=taxa, organisms=sorted(organisms | set(taxa.values())), strains=sorted(strains), accessions=accs, accession_digits=digits)


PROCESSED = [  # (regex on file name, kind) in order of preference
    (r"_FeatureCounts_Unnormalized_Counts_GLbulkRNAseq\.csv$", "counts"),
    (r"_Raw_Gene_Counts\.txt(\.gz)?$", "counts"),
    (r"_FeatureCounts_GLbulkRNAseq\.tsv$", "featurecounts"),
    (r"_array_differential_expression_GLmicroarray\.csv$", "array"),
]


def describe(osd, cache=None):
    """Everything useful about one study: organism(s), taxids, strain, assay, processed tables, genome hints."""
    n = osd_number(osd)
    meta, files = get_json(f"{API}/osd/meta/{n}", cache), get_json(f"{API}/osd/files/{n}", cache)
    info = parse_meta(meta, n)
    flist = files["studies"].get(f"OSD-{n}", {}).get("study_files", [])
    names = {f["file_name"]: f for f in flist}
    tables = []
    for rx, kind in PROCESSED:
        tables += [(nm, kind) for nm in names if re.search(rx, nm)]
    info.update(osd=n, files=len(flist), tables=tables,
                raw_only=not tables, urls={nm: BASE + f["remote_url"] for nm, f in names.items()})
    return info


def table_ids(info, max_bytes=MAX_TABLE_BYTES, entrez=False):
    """Gene identifiers from the study's best processed table, or None.  -> (ids, file_name, kind).

    With `entrez`, array tables give their Entrez gene ids (blank/NA dropped) instead of probe set ids."""
    for name, kind in info["tables"]:
        size = int(_head_size(info["urls"][name]) or 0)
        if size > max_bytes:
            continue
        raw = _get(info["urls"][name])
        if name.endswith(".gz"):
            raw = gzip.decompress(raw)
        delim = "\t" if re.sub(r"\.gz$", "", name).endswith(("tsv", "txt")) else ","
        rows = list(csv.reader(io.StringIO(raw.decode("utf-8-sig", "replace")), delimiter=delim))
        if kind == "array":
            col = rows[0].index("ENTREZID" if entrez else "ProbesetID")
            return [r[col] for r in rows[1:] if len(r) > col and r[col] not in ("", "NA")], name, kind
        return [r[0] for r in rows[1:] if r and r[0]], name, kind
    return None


def _head_size(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "gene-lingo/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.headers.get("Content-Length")
    except urllib.error.URLError:
        return None


# ---- genome choice --------------------------------------------------------------

def candidate_genomes(info, cache=None, index=None):
    """[(acc, row_or_None, why)] best guess first: assembly named in the metadata, then taxon/strain matches."""
    idx = index if index is not None else genomes.assembly_index(cache)
    out, seen = [], set()

    def add(acc, row, why):
        if acc not in seen:
            seen.add(acc)
            out.append((acc, row, why))

    for acc in info["accessions"]:
        add(acc, genomes.lookup(acc, index=idx), "assembly named in study metadata")
    for d in info["accession_digits"]:
        row = genomes.lookup("GCA_" + d, index=idx)
        add(row["acc"] if row else f"GCA_{d}.1", row, "assembly named in study metadata (annotation file name)")
    for name in info["organisms"]:
        for strain in info["strains"] or [None]:
            q = f"{name} {strain}" if strain and strain.lower() not in name.lower() else name if strain is None else None
            for row in (genomes.find_assemblies(query=q, index=idx, limit=2) if q and len(q.split()) > 2 else []):
                add(row["acc"], row, f"matches organism+strain ({q})")
            # strain text that names the parent ("chi3339 (mouse-passaged derivative of SL1344)")
            for tok in re.findall(r"[A-Za-z]*\d[A-Za-z0-9\-]*", strain or ""):
                for row in genomes.find_assemblies(query=f"{' '.join(name.split()[:2])} {tok}", index=idx, limit=1):
                    if genomes._strain_match(row["strain"], [tok.lower()]) == 3:
                        add(row["acc"], row, f"matches organism+strain ({name} {tok})")
    unremarkable = []   # taxon matches that are neither a named strain nor a reference/representative genome
    for taxid in info["taxa"]:
        for row in genomes.find_assemblies(taxid=taxid, index=idx, limit=3):
            if row["category"] in genomes.CATEGORY_RANK or row["taxid"] != taxid:
                add(row["acc"], row, "NCBI taxon of the study" if row["taxid"] == taxid else "reference/representative genome of the species")
            else:
                unremarkable.append(row)
    for name in info["organisms"]:
        for row in genomes.find_assemblies(query=" ".join(name.split()[:2]), index=idx, limit=2):
            add(row["acc"], row, "reference/representative genome of the species")
    for row in unremarkable:   # an arbitrary strain of the taxon is a last resort
        add(row["acc"], row, "NCBI taxon of the study")
    return out


def id_variants(i):
    """Spellings of a table id worth trying as a locus tag: itself, without an Affymetrix '_at' suffix, its first token."""
    base = re.sub(r"_(?:[sxari]_)?at$", "", i)
    return dict.fromkeys([i, base, base.split("_")[0]])


def choose_genome(info, cache=None, work_cache=None, index=None, ids=None, entrez=None):
    """Pick the genome for one study -> dict(acc, label, species, confidence, evidence, matched, n_ids).

    With identifiers (`ids`; `entrez` for arrays annotated with Entrez gene ids) each candidate is scored by how
    many of them it recognises, which separates strains and annotation versions; without, the best
    metadata/taxon match is returned, graded stated / strain / species / taxon.
    """
    cands = candidate_genomes(info, cache, index)
    if not cands:
        return None
    cache = work_cache or cache or genomes.default_cache()
    best, scored = None, []
    if ids or entrez:
        sample = [i for i in (ids or []) if i][:3000]
        esample = [e for e in (entrez or []) if e][:3000]
        for acc, row, why in cands[:6]:
            try:
                gff, faa = bacteria.fetch_genome(acc, cache, {acc: dict(acc=acc)})
            except (urllib.error.URLError, AttributeError, OSError):
                continue
            g = bacteria.Proteome(acc, gff, faa)
            by_id = sum(any(g.resolve(v) for v in id_variants(i)) for i in sample), len(sample)
            by_entrez = sum(any(x in g.geneid for x in e.split("|")) for e in esample), len(esample)
            hit, n = max(by_id, by_entrez, key=lambda hn: hn[0] / hn[1] if hn[1] else 0)
            if n:
                scored.append((hit / n, hit, n, acc, row, why))
        best = max(scored, key=lambda x: x[0], default=None)
        # candidates are in order of how well their *names* match the study; take the first that matches the
        # identifiers about as well as the best one (multi-strain arrays recognise several genomes)
        for frac, hit, n, acc, row, why in scored:
            if best and frac >= max(MIN_MATCH, best[0] - TIE_MARGIN):
                return _result(acc, row, "validated", f"{hit}/{n} identifiers found in genome ({why})", frac,
                               len(ids or entrez), info)
        best = (best[1], best[3], best[4], best[5], best[2]) if best else None
    acc, row, why = cands[0]
    conf = "stated" if why.startswith("assembly named") else "strain" if "strain" in why else "species" if "species" in why else "taxon"
    note = why + (f"; no candidate matched the identifiers (best {best[0]}/{best[4]})" if best else "")
    return _result(acc, row, conf, note, None, len(ids) if ids else 0, info)


def _result(acc, row, confidence, evidence, matched, n_ids, info):
    if row:
        spec = genomes.spec_from_row(row)
    else:   # assembly not in the RefSeq index (GenBank-only): name it from the study
        org = (info["organisms"] or [acc])[0]
        sp = org.split()
        spec = dict(label=" ".join([f"{sp[0][0]}. {sp[1]}" if len(sp) > 1 else org, *info["strains"][:1]]),
                    species=" ".join(sp[:2]), osds=[])
    return dict(acc=acc, label=spec["label"], species=spec["species"], confidence=confidence,
                evidence=evidence, matched=matched, n_ids=n_ids)


def is_bacterial(info, taxa, names):
    """Organism is a bacterial/archaeal taxon (`taxa`: taxids, `names`: lower-case 'genus species') known to the assembly index."""
    return bool(set(info["taxa"]) & taxa) or any(" ".join(o.lower().split()[:2]) in names for o in info["organisms"])


SURVEY_COLS = ["osd", "title", "organism", "taxid", "strain", "technology", "processed_table", "n_ids", "genome",
               "genome_label", "confidence", "id_match", "evidence"]


def survey(cache=None, work_cache=None, progress=None, studies=None):
    """Genome choice for every bacterial/archaeal transcription study in OSDR -> list of row dicts."""
    idx = genomes.assembly_index(cache)
    species_taxa = {r["taxid"] for r in idx} | {r["species_taxid"] for r in idx}
    species_names = {genomes.species_name(r["organism"]).lower() for r in idx}
    note = progress or (lambda m: None)
    rows = []
    done = set()
    for s in (studies if studies is not None else search_studies()):
        if "transcription" not in s["Study Assay Measurement Type"]:
            continue
        osd = osd_number(s["Accession"])
        if osd in done:
            continue
        done.add(osd)
        try:
            info = describe(osd, os.path.join(cache or genomes.default_cache(), "osdr_json"))
        except (OSError, KeyError, IndexError, json.JSONDecodeError):    # OSError covers URLError
            continue
        if not is_bacterial(info, species_taxa, species_names):
            continue
        note(f"OSD-{osd}: {info['title'][:60]}")
        got = None
        try:
            got = table_ids(info)
        except (urllib.error.URLError, ValueError, OSError):
            pass
        ids = got[0] if got else None
        ent = table_ids(info, entrez=True)[0] if got and got[2] == "array" else None
        pick = choose_genome(info, cache, work_cache, idx, ids, ent)
        rows.append(dict(
            osd=f"OSD-{osd}", title=info["title"], organism="; ".join(info["organisms"]), taxid=";".join(info["taxa"]),
            strain="; ".join(info["strains"]), technology=s.get("Study Assay Technology Type", "").split("     ")[0].strip(),
            processed_table=got[1] if got else "", n_ids=len(ids) if ids else 0,
            genome=pick["acc"] if pick else "", genome_label=pick["label"] if pick else "",
            confidence=pick["confidence"] if pick else "none",
            id_match=f"{pick['matched']:.0%}" if pick and pick["matched"] is not None else "",
            evidence=pick["evidence"] if pick else "no genome found"))
    return rows


def write_survey(rows, path):
    with open(path, "w") as fh:
        w = csv.DictWriter(fh, SURVEY_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def genomes_for(osds, cache=None, work_cache=None, progress=None, index=None):
    """Resolve several studies at once -> (registry {acc: spec with osds}, per-study choices {osd: pick}).

    Studies that share an assembly share a registry entry (so their identifiers map to the same genes)."""
    idx = index if index is not None else genomes.assembly_index(cache)
    note = progress or (lambda m: None)
    registry, picks = {}, {}
    for osd in osds:
        n = osd_number(osd)
        note(f"reading OSD-{n}")
        info = describe(n, os.path.join(cache or genomes.default_cache(), "osdr_json"))
        got = None
        try:
            got = table_ids(info)
        except (urllib.error.URLError, ValueError, OSError):
            pass
        ent = table_ids(info, entrez=True)[0] if got and got[2] == "array" else None
        pick = choose_genome(info, cache, work_cache, idx, got[0] if got else None, ent)
        if pick is None:
            raise ValueError(f"OSD-{n}: could not find a genome for {info['organisms'] or 'this study'}")
        pick["info"] = info
        picks[n] = pick
        spec = registry.setdefault(pick["acc"], dict(label=pick["label"], species=pick["species"], acc=pick["acc"], osds=[]))
        spec["osds"].append(n)
    return registry, picks
