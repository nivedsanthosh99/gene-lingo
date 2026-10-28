"""Locus-tag mapping between bacterial genomes (built for the GeneLab/OSDR datasets, works for any genomes).

Orthology is computed from protein sequences (reciprocal best DIAMOND hits), never from names, and every
result is keyed on the locus tags of the annotation the data was processed with.  Everything is expressed
through a hub genome and also directly for every pair of genomes.  Needs the `diamond` binary on PATH.

`GENOMES` is the curated registry for the OSDR studies of five species; every function that takes
`genomes=` accepts any registry {key: {label, species, acc, osds}} (see genomes.py and osdr.py for
building one from names or from OSDR study ids).
"""
import csv
import gzip
import itertools
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.parse
import urllib.request

FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# key -> species, assembly the OSDR processing used (or, for raw-only studies, the
# matching reference strain), OSD studies whose locus tags it covers
GENOMES = {
    "paer_pao1": dict(label="P. aeruginosa PAO1", species="Pseudomonas aeruginosa", acc="GCF_000006765.1", osds=[14, 15]),
    "paer_pa14": dict(label="P. aeruginosa PA14", species="Pseudomonas aeruginosa", acc="GCF_000014625.1", osds=[554]),
    "bsub_168": dict(label="B. subtilis 168", species="Bacillus subtilis", acc="GCF_000009045.1", osds=[138, 185, 215, 216], raw_only=[215, 216]),
    "saur_uams1": dict(label="S. aureus UAMS-1", species="Staphylococcus aureus", acc="GCA_000788115.1", osds=[145]),
    "ecol_k12": dict(label="E. coli K-12 MG1655", species="Escherichia coli", acc="GCF_000005845.2", osds=[43]),
    "ecol_39282": dict(label="E. coli ATCC 4157 (NCTC 86)", species="Escherichia coli", acc="GCF_900699165.1", osds=[95]),
    "ecol_rel606": dict(label="E. coli B REL606", species="Escherichia coli", acc="GCF_000017985.1", osds=[728], raw_only=[728]),
    "rrub_11170": dict(label="R. rubrum ATCC 11170", species="Rhodospirillum rubrum", acc="GCF_000013085.1", osds=[31, 39]),
}
HUB = "bsub_168"


def genome_for_osd(osd):
    osd = int(re.sub(r"\D", "", str(osd)))
    for key, g in GENOMES.items():
        if osd in g["osds"]:
            return key
    raise KeyError(f"OSD-{osd} is not in GENOMES")


# ---- download ----------------------------------------------------------------

def fetch_genome(key, cache, genomes=None):
    """Download genomic.gff.gz and protein.faa.gz for genomes[key] (default GENOMES); return their paths.

    A spec with local `gff` and `faa` paths (gzipped) is used as is, without any download.
    """
    spec = (genomes or GENOMES)[key]
    if spec.get("gff") and spec.get("faa"):
        return spec["gff"], spec["faa"]
    acc = spec["acc"]
    folder = f"{FTP}/{acc[:3]}/{acc[4:7]}/{acc[7:10]}/{acc[10:13]}/"
    with urllib.request.urlopen(folder, timeout=60) as r:
        name = re.search(rf'href="({acc}_[^"/]+)/"', r.read().decode()).group(1)
    os.makedirs(cache, exist_ok=True)
    paths = []
    for ext in ("genomic.gff.gz", "protein.faa.gz"):
        dest = os.path.join(cache, f"{name}_{ext}")
        if not os.path.exists(dest):
            urllib.request.urlretrieve(f"{folder}{name}/{name}_{ext}", dest)
        paths.append(dest)
    return tuple(paths)


# ---- proteomes ---------------------------------------------------------------

def _attrs(s):
    return {k: urllib.parse.unquote(v) for k, v in (kv.split("=", 1) for kv in s.split(";") if "=" in kv)}


def _read_fasta(path):
    seqs, name = {}, None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


class GeneTable:
    """Protein-coding genes of one genome: locus_tag -> {protein_id, symbol, product, old, geneid, ...}."""

    def __init__(self, key, cds):
        self.key, self.cds = key, cds
        self.alias = {}  # alternative locus tags / old tags / punctuation variants -> canonical tag
        for tag, rec in cds.items():
            for alt in (tag, *rec["old"], tag.replace("_", "")):
                self.alias.setdefault(alt, tag)
        self.geneid = {r["geneid"]: t for t, r in cds.items() if r["geneid"]}

    def resolve(self, name):
        """Map a locus tag as written in a dataset (current, old, or underscore-less) to the canonical one."""
        return self.alias.get(name) or self.alias.get(name.replace("_", ""))


class Proteome(GeneTable):
    """GeneTable parsed from an NCBI GFF + protein FASTA; each record also has the protein `seq`."""

    def __init__(self, key, gff_path, faa_path):
        seqs = _read_fasta(faa_path)
        gene_info, cds = {}, {}
        with gzip.open(gff_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] not in ("gene", "pseudogene", "CDS"):
                    continue
                a = _attrs(f[8])
                if f[2] == "CDS":
                    pid, tag = a.get("protein_id"), a.get("locus_tag")
                    if not (pid and tag) or pid not in seqs or tag in cds:
                        continue
                    cds[tag] = dict(protein_id=pid, seq=seqs[pid], product=a.get("product", ""))
                else:
                    gene_info[a.get("locus_tag")] = a
        for tag, rec in cds.items():
            g = gene_info.get(tag, {})
            rec["symbol"] = g.get("gene", "")
            rec["old"] = [x for x in g.get("old_locus_tag", "").split(",") if x]
            m = re.search(r"GeneID:(\d+)", g.get("Dbxref", ""))
            rec["geneid"] = m.group(1) if m else ""
        super().__init__(key, cds)

    def write_fasta(self, path):
        with open(path, "w") as fh:
            for tag, r in self.cds.items():
                fh.write(f">{tag}\n{r['seq']}\n")


# ---- probe/reporter ids for array studies ------------------------------------

def probe_to_locus(osd, probe, proteome):
    """Locus tag for an array probe id in the OSDs that report probes instead of locus tags."""
    osd = int(re.sub(r"\D", "", str(osd)))
    if osd in (14, 15):                       # Affymetrix P. aeruginosa chip, PAO1: 'PA0001_dnaA_at'
        m = re.match(r"(PA\d{4})(?:_|$)", probe)
        return proteome.resolve(m.group(1)) if m else None
    if osd in (31, 39):                       # SCKCEN R. rubrum array: 'RruA0001' -> Rru_A0001
        m = re.match(r"Rru(A\d{4})$", probe)
        return proteome.resolve("Rru_" + m.group(1)) if m else None
    return None


def entrez_to_locus(entrez, proteome):
    """OSD-43 style 'id1|id2' Entrez ids -> locus tags found in this genome."""
    return [proteome.geneid[e] for e in entrez.split("|") if e in proteome.geneid]


# ---- orthology ---------------------------------------------------------------

def _diamond(query_fa, db_prefix, out, threads, sensitivity):
    subprocess.run(["diamond", "blastp", "-q", query_fa, "-d", db_prefix, "-o", out, "--threads", str(threads),
                    "--evalue", "1e-5", "-k", "5", f"--{sensitivity}", "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore"],
                   check=True)


def _best_hits(path, min_cov, min_id):
    """Best hit per query among alignments covering >= min_cov of both proteins."""
    best = {}
    with open(path) as fh:
        for q, s, pid, ln, ql, sl, bs in csv.reader(fh, delimiter="\t"):
            ln, ql, sl, bs, pid = int(ln), int(ql), int(sl), float(bs), float(pid)
            if ln / ql < min_cov or ln / sl < min_cov or pid < min_id:
                continue
            if q not in best or bs > best[q][1]:
                best[q] = (s, bs, pid)
    return best


def rbh(a, b, work, threads=4, min_cov=0.5, min_id=0.0, sensitivity="very-sensitive"):
    """Reciprocal best hits between Proteomes a and b -> {tag_a: (tag_b, pident)}."""
    os.makedirs(work, exist_ok=True)
    files = {}
    for p in (a, b):
        fa = os.path.join(work, f"{p.key}.faa")
        p.write_fasta(fa)
        subprocess.run(["diamond", "makedb", "--in", fa, "-d", os.path.join(work, p.key), "--quiet"], check=True)
        files[p.key] = fa
    ab, ba = (os.path.join(work, f"{x.key}_vs_{y.key}.tsv") for x, y in ((a, b), (b, a)))
    _diamond(files[a.key], os.path.join(work, b.key), ab, threads, sensitivity)
    _diamond(files[b.key], os.path.join(work, a.key), ba, threads, sensitivity)
    fwd, rev = _best_hits(ab, min_cov, min_id), _best_hits(ba, min_cov, min_id)
    return {q: (s, pid) for q, (s, _, pid) in fwd.items() if rev.get(s, (None,))[0] == q}


def _cached_rbh(cache, specs, a, b, work, **kw):
    """rbh() with results kept in <cache>/rbh keyed by assembly accession and thresholds, so adding a
    genome to a set only computes the new pairs.  Local-file genomes are not cached."""
    accs = [specs[k].get("acc") for k in (a.key, b.key)]
    if not all(accs) or any(specs[k].get("gff") for k in (a.key, b.key)):
        return rbh(a, b, work, **kw)
    sig = f"c{kw.get('min_cov', 0.5)}_i{kw.get('min_id', 0.0)}_{kw.get('sensitivity', 'very-sensitive')}"
    flip = accs[0] > accs[1]
    x, y = (b, a) if flip else (a, b)
    path = os.path.join(cache, "rbh", f"{min(accs)}__{max(accs)}__{sig}.tsv")
    if os.path.exists(path):
        with open(path) as fh:
            rows = [(r[0], r[1], float(r[2])) for r in csv.reader(fh, delimiter="\t")]
    else:
        res = rbh(x, y, work, **kw)
        rows = [(t, v[0], v[1]) for t, v in res.items()]
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    return {t: (u, pid) for t, u, pid in rows} if not flip else {u: (t, pid) for t, u, pid in rows}


def build(cache, work, keys=None, hub=HUB, genomes=None, hub_only=False, progress=None, **kw):
    """Run everything: returns (proteomes, hub_map, pair_map).

    hub_map[key][hub_tag] = (tag_in_key, pident); pair_map[(a, b)][tag_a] = (tag_b, pident)
    `genomes` is a registry like GENOMES; with hub_only just the hub-vs-genome pairs are computed
    (n-1 instead of n(n-1)/2 comparisons).  `progress(msg)` is called as work proceeds.
    """
    if not shutil.which("diamond"):
        raise RuntimeError("the `diamond` binary is required (apt install diamond-aligner)")
    specs = genomes or GENOMES
    keys = list(keys or specs)
    note = progress or (lambda msg: None)
    prot = {}
    for k in keys:
        note(f"loading {specs[k].get('label', k)}")
        prot[k] = Proteome(k, *fetch_genome(k, cache, specs))
    pair = {}
    todo = [(hub, k) for k in keys if k != hub] if hub_only else list(itertools.combinations(keys, 2))
    for i, (a, b) in enumerate(todo, 1):
        note(f"orthologs {i}/{len(todo)}: {specs[a].get('label', a)} vs {specs[b].get('label', b)}")
        pair[(a, b)] = _cached_rbh(cache, specs, prot[a], prot[b], os.path.join(work, f"{a}__{b}"), **kw)
    hub_map = {}
    for k in keys:
        if k == hub:
            continue
        if (hub, k) in pair:
            hub_map[k] = pair[(hub, k)]
        else:
            hub_map[k] = {v[0]: (h, v[1]) for h, v in pair[(k, hub)].items()}
    return prot, hub_map, pair


def build_set(registry, outdir, cache, hub=None, hub_only=False, threads=4, progress=None, **kw):
    """build() for any genome registry, then write_tables() into `outdir`.  The hub defaults to the first genome."""
    if len(registry) < 2:
        raise ValueError("need at least two genomes")
    hub = hub if hub in registry else next(iter(registry))
    with tempfile.TemporaryDirectory() as work:
        prot, hub_map, pair = build(cache, work, hub=hub, genomes=registry, hub_only=hub_only,
                                    progress=progress, threads=threads, **kw)
    write_tables(outdir, prot, hub_map, pair, hub, registry)
    return prot, hub_map, pair


def hub_rows(prot, hub_map, hub=HUB):
    """One row per hub CDS: hub tag, symbol, product, then one ortholog tag per other genome."""
    keys = list(hub_map)
    for tag, r in prot[hub].cds.items():
        yield [tag, r["symbol"], r["product"]] + [hub_map[k].get(tag, ("",))[0] for k in keys]


def core(hub_map, species=None, genomes=None):
    """Hub tags having an ortholog in every genome (or, with `species`, every listed species' genomes)."""
    keys = [k for k in hub_map if species is None or (genomes or GENOMES)[k]["species"] in species]
    sets = [set(hub_map[k]) for k in keys]
    return sorted(set.intersection(*sets)) if sets else []


def species_core(hub_map, genomes=None):
    """Hub tags that have an ortholog in at least one genome of each species (strain-tolerant core)."""
    specs = genomes or GENOMES
    by_species = {}
    for k in hub_map:
        by_species.setdefault(specs[k]["species"], set()).update(hub_map[k])
    return sorted(set.intersection(*by_species.values())) if by_species else []


def write_tables(outdir, prot, hub_map, pair, hub=HUB, genomes=None):
    """hub.tsv (one row per hub CDS), core.tsv, genes.tsv.gz, genomes.json and pairs/<a>__<b>.tsv."""
    specs = genomes or GENOMES
    os.makedirs(os.path.join(outdir, "pairs"), exist_ok=True)
    write_genes(outdir, prot)
    with open(os.path.join(outdir, "genomes.json"), "w") as fh:
        json.dump({"hub": hub, "genomes": {k: specs[k] for k in prot}}, fh, indent=1)
    keys = list(hub_map)
    with open(os.path.join(outdir, "hub.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"{hub}_locus_tag", "symbol", "product", *keys])
        w.writerows(hub_rows(prot, hub_map, hub))
    with open(os.path.join(outdir, "core.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"{hub}_locus_tag", "symbol", "product", *keys])
        in_core = set(species_core(hub_map, specs))
        w.writerows(r for r in hub_rows(prot, hub_map, hub) if r[0] in in_core)
    for (a, b_), m in pair.items():
        with open(os.path.join(outdir, "pairs", f"{a}__{b_}.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([a, b_, "pident"])
            w.writerows((x, y, f"{pid:.1f}") for x, (y, pid) in m.items())


def write_genes(outdir, prot):
    """genes.tsv.gz: every protein-coding gene of every genome, with its alternative names."""
    with gzip.open(os.path.join(outdir, "genes.tsv.gz"), "wt") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "locus_tag", "old_locus_tags", "symbol", "product", "protein_id", "geneid"])
        for key, p in prot.items():
            for tag, r in p.cds.items():
                w.writerow([key, tag, ";".join(r["old"]), r["symbol"], r["product"], r["protein_id"], r["geneid"]])


def load_genes(outdir):
    """Read genes.tsv.gz -> {genome key: GeneTable} (no sequences, no network needed)."""
    cds = {}
    with gzip.open(os.path.join(outdir, "genes.tsv.gz"), "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cds.setdefault(r["genome"], {})[r["locus_tag"]] = dict(
                old=[x for x in r["old_locus_tags"].split(";") if x], symbol=r["symbol"], product=r["product"],
                protein_id=r["protein_id"], geneid=r["geneid"])
    return {k: GeneTable(k, v) for k, v in cds.items()}


def load_pair_scored(outdir, a, b_):
    """Read pairs/<a>__<b>.tsv (either order) -> {tag_a: (tag_b, pident)}."""
    for x, y, swap in ((a, b_, False), (b_, a, True)):
        path = os.path.join(outdir, "pairs", f"{x}__{y}.tsv")
        if os.path.exists(path):
            with open(path) as fh:
                rows = list(csv.reader(fh, delimiter="\t"))[1:]
            return {r[1]: (r[0], float(r[2])) for r in rows} if swap else {r[0]: (r[1], float(r[2])) for r in rows}
    raise KeyError(f"no mapping between {a} and {b_} in {outdir}")


def load_pair(outdir, a, b_):
    """Read pairs/<a>__<b>.tsv (either order) -> {tag_a: tag_b}."""
    return {k: v[0] for k, v in load_pair_scored(outdir, a, b_).items()}


def resolve_input(name, key, genes, genomes=None):
    """Canonical locus tags (usually one) for a dataset identifier of genome `key`.

    Accepts current/old locus tags, array probe ids of the OSDs on that genome, and Entrez gene ids
    (the 'id1|id2' form in the OSD-43 tables).
    """
    name = name.strip()
    g = genes[key]
    if hit := g.resolve(name):
        return [hit]
    for osd in (genomes or GENOMES).get(key, {}).get("osds", []):
        if hit := probe_to_locus(osd, name, g):
            return [hit]
    return entrez_to_locus(name, g)


def detect_genome(names, genes, candidates=None, genomes=None):
    """Genome key in which the most of `names` resolve (None if nothing resolves anywhere)."""
    names = [n for n in names if n.strip()][:500]
    scores = {k: sum(bool(resolve_input(n, k, genes, genomes)) for n in names) for k in (candidates or genes)}
    best = max(scores, key=scores.get)
    return best if scores[best] else None


class Mapper:
    """Translate dataset identifiers between genomes using the tables written by `build`.

    Reads genomes.json for the genome registry and hub (falls back to GENOMES / HUB for older tables).
    Genome arguments may be a key, a label, or an OSD study ("OSD-554", 554).
    """

    def __init__(self, outdir):
        self.outdir = outdir
        self.genes = load_genes(outdir)
        path = os.path.join(outdir, "genomes.json")
        if os.path.exists(path):
            with open(path) as fh:
                meta = json.load(fh)
            self.genomes, self.hub = meta["genomes"], meta["hub"]
        else:
            self.genomes, self.hub = {k: GENOMES[k] for k in self.genes}, HUB
        self._pairs, self._hub_of = {}, {}

    def key_for(self, x):
        """Genome key for a key, label or OSD study id."""
        if x in self.genomes:
            return x
        for k, g in self.genomes.items():
            if x == g.get("label"):
                return k
        m = re.fullmatch(r"(?:OSD-?|GLDS-?)?(\d+)", str(x).strip(), re.I)
        if m:
            for k, g in self.genomes.items():
                if int(m.group(1)) in g.get("osds", []):
                    return k
        raise KeyError(f"unknown genome or study: {x}")

    def pair(self, a, b_):
        """{tag_a: (tag_b, pident)}; pairs that were not computed (hub-only builds) go through the hub."""
        if (a, b_) not in self._pairs:
            try:
                self._pairs[(a, b_)] = load_pair_scored(self.outdir, a, b_)
            except KeyError:
                if self.hub in (a, b_):
                    raise
                to_b = load_pair_scored(self.outdir, self.hub, b_)          # hub tag -> (tag in b, pident)
                self._pairs[(a, b_)] = {t: (to_b[h][0], min(p, to_b[h][1]))
                                        for t, (h, p) in load_pair_scored(self.outdir, a, self.hub).items() if h in to_b}
        return self._pairs[(a, b_)]

    def hub_tag(self, key, tag):
        if key == self.hub:
            return tag
        if key not in self._hub_of:
            self._hub_of[key] = {t: v[0] for t, v in self.pair(key, self.hub).items()}
        return self._hub_of[key].get(tag, "")

    def resolve(self, name, key):
        return resolve_input(name, key, self.genes, self.genomes)

    def detect(self, names, candidates=None):
        return detect_genome(names, self.genes, candidates, self.genomes)

    def map(self, names, src, dsts):
        """One row per (input, resolved gene) with the hub gene and, per target genome, tag + % identity."""
        src, dsts = self.key_for(src), [self.key_for(d) for d in dsts]
        rows = []
        for name in names:
            tags = self.resolve(name, src) or [None]
            for tag in tags:
                rec = self.genes[src].cds.get(tag, {})
                row = {"input": name, "locus_tag": tag or "", "symbol": rec.get("symbol", ""),
                       "product": rec.get("product", ""), "hub_gene": self.hub_tag(src, tag) if tag else ""}
                for d in dsts:
                    hit = None if not tag else (tag, 100.0) if d == src else self.pair(src, d).get(tag)
                    row[d] = hit[0] if hit else ""
                    row[f"{d}_pident"] = round(hit[1], 1) if hit else None
                    row[f"{d}_symbol"] = self.genes[d].cds[hit[0]]["symbol"] if hit else ""
                rows.append(row)
        return rows

    def hub_table(self):
        """The hub.tsv rows as dicts (one per hub gene); the first column is called `locus_tag`."""
        with open(os.path.join(self.outdir, "hub.tsv")) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        first = next(iter(rows[0])) if rows else None
        return [{("locus_tag" if k == first else k): v for k, v in r.items()} for r in rows]


def translate(tags, src, dst, pair_map, proteome=None):
    """Map `src` locus tags to `dst` -> {tag: tag_or_None} using a {tag_src: tag_dst} dict.

    `proteome` (of src) additionally resolves old/underscore-less spellings.
    """
    if src == dst:
        return {t: (proteome.resolve(t) if proteome else t) for t in tags}
    out = {}
    for t in tags:
        c = proteome.resolve(t) if proteome else t
        out[t] = pair_map.get(c) if c else None
    return out
