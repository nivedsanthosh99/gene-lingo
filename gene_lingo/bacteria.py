"""Locus-tag mapping between the bacterial genomes used by GeneLab/OSDR datasets.

Orthology is computed from protein sequences (reciprocal best DIAMOND hits), never
from names, and every result is keyed on the locus tags that appear in the OSDR
count tables.  Everything is expressed through a hub genome (B. subtilis 168) and
also directly for every pair of genomes.  Needs the `diamond` binary on PATH.
"""
import csv
import gzip
import itertools
import os
import re
import shutil
import subprocess
import urllib.parse
import urllib.request

FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# key -> species, assembly the OSDR processing used (or, for raw-only studies, the
# matching reference strain), OSD studies whose locus tags it covers
GENOMES = {
    "paer_pao1": dict(species="Pseudomonas aeruginosa", acc="GCF_000006765.1", osds=[14, 15]),
    "paer_pa14": dict(species="Pseudomonas aeruginosa", acc="GCF_000014625.1", osds=[554]),
    "bsub_168": dict(species="Bacillus subtilis", acc="GCF_000009045.1", osds=[138, 185, 215, 216]),
    "saur_uams1": dict(species="Staphylococcus aureus", acc="GCA_000788115.1", osds=[145]),
    "ecol_k12": dict(species="Escherichia coli", acc="GCF_000005845.2", osds=[43]),
    "ecol_39282": dict(species="Escherichia coli", acc="GCF_900699165.1", osds=[95]),
    "ecol_rel606": dict(species="Escherichia coli", acc="GCF_000017985.1", osds=[728]),
    "rrub_11170": dict(species="Rhodospirillum rubrum", acc="GCF_000013085.1", osds=[31, 39]),
}
HUB = "bsub_168"


def genome_for_osd(osd):
    osd = int(re.sub(r"\D", "", str(osd)))
    for key, g in GENOMES.items():
        if osd in g["osds"]:
            return key
    raise KeyError(f"OSD-{osd} is not in GENOMES")


# ---- download ----------------------------------------------------------------

def fetch_genome(key, cache):
    """Download genomic.gff.gz and protein.faa.gz for GENOMES[key]; return their paths."""
    acc = GENOMES[key]["acc"]
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


class Proteome:
    """CDS records of one genome: locus_tag -> {protein_id, seq, symbol, product, old, geneid}."""

    def __init__(self, key, gff_path, faa_path):
        self.key = key
        seqs = _read_fasta(faa_path)
        gene_info = {}
        self.cds = {}
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
                    if not (pid and tag) or pid not in seqs or tag in self.cds:
                        continue
                    self.cds[tag] = dict(protein_id=pid, seq=seqs[pid], product=a.get("product", ""))
                else:
                    gene_info[a.get("locus_tag")] = a
        self.alias = {}  # alternative locus tags / old tags / punctuation variants -> canonical tag
        for tag, rec in self.cds.items():
            g = gene_info.get(tag, {})
            rec["symbol"] = g.get("gene", "")
            rec["old"] = [x for x in g.get("old_locus_tag", "").split(",") if x]
            m = re.search(r"GeneID:(\d+)", g.get("Dbxref", ""))
            rec["geneid"] = m.group(1) if m else ""
            for alt in (tag, *rec["old"], tag.replace("_", "")):
                self.alias.setdefault(alt, tag)
        self.geneid = {r["geneid"]: t for t, r in self.cds.items() if r["geneid"]}

    def resolve(self, name):
        """Map a locus tag as written in a dataset (current, old, or underscore-less) to the canonical one."""
        return self.alias.get(name) or self.alias.get(name.replace("_", ""))

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


def build(cache, work, keys=None, hub=HUB, **kw):
    """Run everything: returns (proteomes, hub_map, pair_map).

    hub_map[key][hub_tag] = (tag_in_key, pident); pair_map[(a, b)][tag_a] = (tag_b, pident)
    """
    if not shutil.which("diamond"):
        raise RuntimeError("the `diamond` binary is required (apt install diamond-aligner)")
    keys = list(keys or GENOMES)
    prot = {k: Proteome(k, *fetch_genome(k, cache)) for k in keys}
    pair = {}
    for a, b in itertools.combinations(keys, 2):
        pair[(a, b)] = rbh(prot[a], prot[b], os.path.join(work, f"{a}__{b}"), **kw)
    hub_map = {}
    for k in keys:
        if k == hub:
            continue
        if (hub, k) in pair:
            hub_map[k] = pair[(hub, k)]
        else:
            hub_map[k] = {v[0]: (h, v[1]) for h, v in pair[(k, hub)].items()}
    return prot, hub_map, pair


def hub_rows(prot, hub_map, hub=HUB):
    """One row per hub CDS: hub tag, symbol, product, then one ortholog tag per other genome."""
    keys = list(hub_map)
    for tag, r in prot[hub].cds.items():
        yield [tag, r["symbol"], r["product"]] + [hub_map[k].get(tag, ("",))[0] for k in keys]


def core(hub_map, species=None):
    """Hub tags having an ortholog in every genome (or, with `species`, every listed species' genomes)."""
    keys = [k for k in hub_map if species is None or GENOMES[k]["species"] in species]
    sets = [set(hub_map[k]) for k in keys]
    return sorted(set.intersection(*sets)) if sets else []


def species_core(hub_map):
    """Hub tags that have an ortholog in at least one genome of each species (strain-tolerant core)."""
    by_species = {}
    for k in hub_map:
        by_species.setdefault(GENOMES[k]["species"], set()).update(hub_map[k])
    return sorted(set.intersection(*by_species.values()))


def write_tables(outdir, prot, hub_map, pair, hub=HUB):
    """hub.tsv (one row per hub CDS), core.tsv, and pairs/<a>__<b>.tsv for every genome pair."""
    os.makedirs(os.path.join(outdir, "pairs"), exist_ok=True)
    keys = list(hub_map)
    with open(os.path.join(outdir, "hub.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"{hub}_locus_tag", "symbol", "product", *keys])
        w.writerows(hub_rows(prot, hub_map, hub))
    with open(os.path.join(outdir, "core.tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([f"{hub}_locus_tag", "symbol", "product", *keys])
        in_core = set(species_core(hub_map))
        w.writerows(r for r in hub_rows(prot, hub_map, hub) if r[0] in in_core)
    for (a, b_), m in pair.items():
        with open(os.path.join(outdir, "pairs", f"{a}__{b_}.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([a, b_, "pident"])
            w.writerows((x, y, f"{pid:.1f}") for x, (y, pid) in m.items())


def load_pair(outdir, a, b_):
    """Read pairs/<a>__<b>.tsv (either order) -> {tag_a: tag_b}."""
    for x, y, swap in ((a, b_, False), (b_, a, True)):
        path = os.path.join(outdir, "pairs", f"{x}__{y}.tsv")
        if os.path.exists(path):
            with open(path) as fh:
                rows = list(csv.reader(fh, delimiter="\t"))[1:]
            return {r[1]: r[0] for r in rows} if swap else {r[0]: r[1] for r in rows}
    raise KeyError(f"no mapping between {a} and {b_} in {outdir}")


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
