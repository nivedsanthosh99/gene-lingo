import gzip
import random
import shutil

import pytest

from gene_lingo import bacteria

AA = "ACDEFGHIKLMNPQRSTVWY"


def prot_seq(rng, n=150):
    return "".join(rng.choice(AA) for _ in range(n))


def mutate(rng, s, frac):
    return "".join(rng.choice(AA) if rng.random() < frac else c for c in s)


def make_genome(tmp_path, key, genes):
    """genes: [(locus_tag, old_tag, symbol, seq)] -> Proteome built from a synthetic GFF + FASTA."""
    gff, faa = tmp_path / f"{key}.gff.gz", tmp_path / f"{key}.faa.gz"
    with gzip.open(gff, "wt") as g, gzip.open(faa, "wt") as f:
        g.write("##gff-version 3\n")
        for i, (tag, old, sym, seq) in enumerate(genes):
            attrs = f"ID=gene-{tag};Name={sym};gene={sym};gbkey=Gene;locus_tag={tag}"
            if old:
                attrs += f";old_locus_tag={old}"
            attrs += f";Dbxref=GeneID:{1000 + i}"
            g.write(f"chr\tRefSeq\tgene\t{i*1000+1}\t{i*1000+450}\t.\t+\t.\t{attrs}\n")
            g.write(f"chr\tRefSeq\tCDS\t{i*1000+1}\t{i*1000+450}\t.\t+\t0\tID=cds-{tag};Parent=gene-{tag};"
                    f"locus_tag={tag};protein_id=P_{tag};product=prod%20{sym}\n")
            f.write(f">P_{tag} prod {sym}\n{seq}\n")
        g.write(f"chr\tRefSeq\tpseudogene\t90000\t90300\t.\t+\t.\tID=gene-PSEUDO;locus_tag=PSEUDO\n")
    return bacteria.Proteome(key, gff, faa)


@pytest.fixture
def genomes(tmp_path):
    rng = random.Random(7)
    base = [prot_seq(rng) for _ in range(4)]
    a = make_genome(tmp_path, "A", [(f"A_RS{i}", f"A_{i}", f"g{i}", s) for i, s in enumerate(base)])
    # B: gene 0 identical, gene 1 diverged, gene 2 missing, an unrelated gene; different tag scheme
    b = make_genome(tmp_path, "B", [("B_001", "", "x0", base[0]), ("B_002", "", "x1", mutate(rng, base[1], 0.3)),
                                   ("B_003", "", "x3", base[3]), ("B_004", "", "zz", prot_seq(rng))])
    return a, b


def test_proteome_parsing_and_resolve(genomes):
    a, _ = genomes
    assert list(a.cds) == ["A_RS0", "A_RS1", "A_RS2", "A_RS3"]      # pseudogene has no CDS -> skipped
    assert a.cds["A_RS1"]["product"] == "prod g1" and a.cds["A_RS1"]["symbol"] == "g1"
    # the old/new locus-tag mismatch that breaks name-based lookups
    assert a.resolve("A_1") == "A_RS1" == a.resolve("A_RS1") == a.resolve("ARS1")
    assert a.resolve("nope") is None
    assert a.geneid["1002"] == "A_RS2"


def test_probe_ids(genomes):
    a, _ = genomes
    assert bacteria.entrez_to_locus("1001|9999|1003", a) == ["A_RS1", "A_RS3"]
    assert bacteria.probe_to_locus(14, "PA0001_dnaA_at", a) is None   # tag not in this toy genome
    pa = bacteria.Proteome.__new__(bacteria.Proteome)
    pa.alias = {"PA0001": "PA0001", "PA0005": "PA0005"}
    assert bacteria.probe_to_locus(14, "PA0001_dnaA_at", pa) == "PA0001"
    assert bacteria.probe_to_locus(15, "PA0005_s_at", pa) == "PA0005"
    assert bacteria.probe_to_locus(15, "ig_901934_901046_at", pa) is None
    assert bacteria.probe_to_locus(554, "PA0001", pa) is None
    rr = bacteria.Proteome.__new__(bacteria.Proteome)
    rr.alias = {"Rru_A0001": "RRU_RS00005"}
    assert bacteria.probe_to_locus(31, "RruA0001", rr) == "RRU_RS00005"


def test_genome_for_osd():
    assert bacteria.genome_for_osd("OSD-554") == "paer_pa14"
    assert bacteria.genome_for_osd(185) == "bsub_168"
    with pytest.raises(KeyError):
        bacteria.genome_for_osd(1)


@pytest.mark.skipif(not shutil.which("diamond"), reason="diamond not installed")
def test_rbh_and_tables(genomes, tmp_path, monkeypatch):
    a, b = genomes
    monkeypatch.setitem(bacteria.GENOMES, "B", dict(species="Toy b", acc="", osds=[]))
    m = bacteria.rbh(a, b, str(tmp_path / "w"), threads=1, sensitivity="sensitive")
    assert {k: v[0] for k, v in m.items()} == {"A_RS0": "B_001", "A_RS1": "B_002", "A_RS3": "B_003"}
    assert m["A_RS0"][1] == 100.0

    out = tmp_path / "out"
    bacteria.write_tables(str(out), {"A": a, "B": b}, {"B": {k: (v[0], v[1]) for k, v in m.items()}}, {("A", "B"): m}, hub="A")
    assert bacteria.load_pair(str(out), "A", "B")["A_RS1"] == "B_002"
    assert bacteria.load_pair(str(out), "B", "A")["B_002"] == "A_RS1"      # reverse direction
    rows = (out / "hub.tsv").read_text().splitlines()
    assert rows[0].split("\t")[-1] == "B" and rows[3].split("\t")[0] == "A_RS2" and rows[3].endswith("\t")
    assert [r.split("\t")[0] for r in (out / "core.tsv").read_text().splitlines()[1:]] == ["A_RS0", "A_RS1", "A_RS3"]

    t = bacteria.translate(["A_1", "A_RS3", "A_RS2", "zzz"], "A", "B", bacteria.load_pair(str(out), "A", "B"), a)
    assert t == {"A_1": "B_002", "A_RS3": "B_003", "A_RS2": None, "zzz": None}


def test_core_helpers():
    hm = {"paer_pao1": {"h1": ("p1", 50), "h2": ("p2", 50)}, "paer_pa14": {"h1": ("q1", 50)},
          "ecol_k12": {"h1": ("e1", 50), "h2": ("e2", 50), "h3": ("e3", 50)}}
    assert bacteria.core(hm) == ["h1"]
    assert bacteria.species_core(hm) == ["h1", "h2"]          # h2 is in at least one genome of each species
