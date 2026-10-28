"""Offline tests for genomes.py (any genome), osdr.py (any OSDR study) and the Mapper's generalisations."""
import gzip
import json
import random
import shutil

import pytest
from helpers import make_genome, mutate, prot_seq

from gene_lingo import bacteria, genomes, osdr


def row(acc, organism, strain="", category="na", level="Complete Genome", taxid="1", species_taxid="1", date="2020-01-01"):
    return dict(acc=acc, category=category, taxid=taxid, species_taxid=species_taxid, organism=organism, strain=strain,
                level=level, asm_name="x", date=date, group="bacteria")


INDEX = [
    row("GCF_000011805.1", "Aliivibrio fischeri ES114", "ES114", "reference genome", taxid="312309", species_taxid="668"),
    row("GCF_000000001.1", "Aliivibrio fischeri", "CG101", taxid="668", species_taxid="668"),
    row("GCF_000013425.1", "Staphylococcus aureus subsp. aureus NCTC 8325", "NCTC 8325", "reference genome", taxid="93061", species_taxid="1280"),
    row("GCF_000000002.1", "Staphylococcus aureus", "", taxid="1280", species_taxid="1280"),
    row("GCF_000017985.1", "Escherichia coli B str. REL606", "REL606", taxid="413997", species_taxid="562"),
]


def test_find_assemblies_ranking():
    q = lambda *a, **k: [r["acc"] for r in genomes.find_assemblies(*a, index=INDEX, **k)]  # noqa: E731
    assert q("Vibrio fischeri ES114")[0] == "GCF_000011805.1"            # genus synonym + strain word
    assert q("Aliivibrio fischeri")[0] == "GCF_000011805.1"              # reference genome outranks other strains
    assert q(taxid=1280)[0] == "GCF_000013425.1"                         # species taxon -> reference, not a strain-less entry
    assert q(taxid=413997) == ["GCF_000017985.1"]                        # strain-level taxon is exact
    assert q("Escherichia coli strain REL606") == ["GCF_000017985.1"]   # the word 'strain' is ignored
    assert q("Nothing here") == []


def test_ranking_regressions():
    idx = INDEX + [
        row("GCF_000007465.2", "Streptococcus mutans UA159", "UA159", "reference genome", taxid="210007", species_taxid="1309"),
        row("GCF_040733985.1", "Streptococcus mutans UA159", "UA159 ykuR deletion", taxid="210007", species_taxid="1309", date="2024-01-01"),
        row("GCF_000000003.1", "Salmonella enterica subsp. enterica", "245-D1", taxid="59201", species_taxid="28901"),
        row("GCF_000006945.2", "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2", "LT2", "reference genome", taxid="99287", species_taxid="28901"),
    ]
    top = lambda **k: genomes.find_assemblies(index=idx, **k)[0]["acc"]  # noqa: E731
    assert top(query="Streptococcus mutans UA159") == "GCF_000007465.2"      # whole-strain match beats a mutant that contains the word
    assert top(taxid=28901) == "GCF_000006945.2"                              # species taxon -> reference genome, not an arbitrary strain


def test_lookup_and_species_label():
    assert genomes.lookup("GCA_000011805.1", index=INDEX)["acc"] == "GCF_000011805.1"   # GCA/GCF share digits
    r = INDEX[2]
    assert genomes.species_name(r["organism"]) == "Staphylococcus aureus"
    assert genomes.label_for(r) == "S. aureus NCTC 8325"


def test_parse_genome_arg(tmp_path):
    key, spec = genomes.parse_genome_arg("GCF_000017985.1", index=INDEX)
    assert key == "GCF_000017985.1" and spec["label"] == "E. coli REL606" and spec["species"] == "Escherichia coli"
    assert genomes.parse_genome_arg("GCA_999999999.1", index=INDEX)[1]["label"] == "GCA_999999999.1"  # unknown accession still usable
    assert genomes.parse_genome_arg("taxid:93061", index=INDEX)[0] == "GCF_000013425.1"
    assert genomes.parse_genome_arg("fischeri es114", index=INDEX)[0] == "GCF_000011805.1"
    k, local = genomes.parse_genome_arg("a.gff.gz,b.faa.gz,My strain", index=INDEX)
    assert k == "My_strain" and local["gff"] == "a.gff.gz" and local["faa"] == "b.faa.gz"
    with pytest.raises(ValueError):
        genomes.parse_genome_arg("no such bug", index=INDEX)


def meta(extra=None):
    tax = lambda n: {"annotationValue": "Vibrio fischeri strain ES114", "termAccession": n}  # noqa: E731
    return {"study": {"OSD-9": {"studies": [{
        "title": "T", "characteristicCategories": [{"@id": "c1", "characteristicType": {"annotationValue": "Organism"}},
                                                   {"@id": "c2", "characteristicType": {"annotationValue": "Strain"}}],
        "materials": {"sources": [{"characteristics": [
            {"category": {"@id": "c1"}, "value": tax("http://purl.obolibrary.org/obo/NCBITaxon_312309")},
            {"category": {"@id": "c2"}, "value": {"annotationValue": "ES114"}},
            {"category": {"@id": "c2"}, "value": {"annotationValue": "Not Applicable"}}]}]},
        "protocols": [{"description": "mapped to GCF_000011805.1_ASM1180v1.fna with gtf Aliivibrio_gca_000009045.ASM904v1.59.gtf",
                       **(extra or {})}]}]}}}


def test_parse_meta():
    m = osdr.parse_meta(meta(), 9)
    assert m["taxa"] == {"312309": "Vibrio fischeri strain ES114"} and m["strains"] == ["ES114"]
    assert m["organisms"] == ["Vibrio fischeri strain ES114"]
    assert m["accessions"] == ["GCF_000011805.1"] and m["accession_digits"] == ["000009045"]
    other = meta()
    other["study"]["OSD-9"]["studies"][0]["materials"]["sources"][0]["characteristics"][0]["value"]["termAccession"] = \
        "http://purl.bioontology.org/ontology/NCBITAXON/287"
    assert osdr.parse_meta(other, 9)["taxa"] == {"287": "Vibrio fischeri strain ES114"}


def test_describe_finds_processed_tables(monkeypatch):
    names = ["GLDS-9_rna_seq_FeatureCounts_Unnormalized_Counts_GLbulkRNAseq.csv", "GLDS-9_rna_seq_processed_md5sum_GLbulkRNAseq.tsv",
             "GLDS-9_rna-seq_GSE1_Raw_Gene_Counts.txt.gz", "GLDS-9_rna_seq_FeatureCounts_rRNArm_GLbulkRNAseq.tsv", "x_R1_raw.fastq.gz"]
    files = {"studies": {"OSD-9": {"study_files": [{"file_name": n, "remote_url": "/d/" + n} for n in names]}}}
    monkeypatch.setattr(osdr, "get_json", lambda url, cache=None: meta() if "/meta/" in url else files)
    info = osdr.describe("OSD-9")
    assert [t[0] for t in info["tables"]] == names[:3][:1] + [names[2]]      # md5sum and rRNArm tsv are not gene tables
    assert not info["raw_only"]
    monkeypatch.setattr(osdr, "get_json", lambda url, cache=None: meta() if "/meta/" in url else {"studies": {"OSD-9": {"study_files": []}}})
    assert osdr.describe(9)["raw_only"]


def test_table_ids(monkeypatch):
    body = gzip.compress(b"Locus Tag\ts1\ts2\nM495_RS00005\t1\t2\nM495_RS00010\t3\t4\n")
    info = {"tables": [("t_Raw_Gene_Counts.txt.gz", "counts")], "urls": {"t_Raw_Gene_Counts.txt.gz": "u"}}
    monkeypatch.setattr(osdr, "_get", lambda url: body)
    monkeypatch.setattr(osdr, "_head_size", lambda url: "10")
    assert osdr.table_ids(info) == (["M495_RS00005", "M495_RS00010"], "t_Raw_Gene_Counts.txt.gz", "counts")
    arr = b'"ENSEMBL","ProbesetID","x"\n"NA","PA0001_dnaA_at","1"\n"NA","ig_1_2_at","2"\n'
    info = {"tables": [("a_array_differential_expression_GLmicroarray.csv", "array")], "urls": {"a_array_differential_expression_GLmicroarray.csv": "u"}}
    monkeypatch.setattr(osdr, "_get", lambda url: arr)
    assert osdr.table_ids(info)[0] == ["PA0001_dnaA_at", "ig_1_2_at"]
    monkeypatch.setattr(osdr, "_head_size", lambda url: str(10 ** 9))          # too big -> skipped, not downloaded
    assert osdr.table_ids(info) is None
    assert list(osdr.id_variants("PA0001_dnaA_at")) == ["PA0001_dnaA_at", "PA0001_dnaA", "PA0001"]


def test_choose_genome_validates_against_identifiers(tmp_path, monkeypatch):
    rng = random.Random(3)
    seqs = [prot_seq(rng) for _ in range(5)]
    right = make_genome(tmp_path, "right", [(f"R_RS{i}", "", f"g{i}", s) for i, s in enumerate(seqs)])
    wrong = make_genome(tmp_path, "wrong", [(f"W_RS{i}", "", f"g{i}", s) for i, s in enumerate(seqs)])
    files = {"GCF_1.1": (tmp_path / "right.gff.gz", tmp_path / "right.faa.gz"), "GCF_2.1": (tmp_path / "wrong.gff.gz", tmp_path / "wrong.faa.gz")}
    monkeypatch.setattr(bacteria, "fetch_genome", lambda key, cache, g=None: files[key])
    rows = [row("GCF_2.1", "Foo bar", "", "reference genome", taxid="10", species_taxid="10"),
            row("GCF_1.1", "Foo bar", "S1", "na", taxid="11", species_taxid="10")]
    info = dict(accessions=[], accession_digits=[], taxa={"10": "Foo bar"}, organisms=["Foo bar"], strains=[])
    ids = [f"R_RS{i}" for i in range(5)]
    got = osdr.choose_genome(info, index=rows, ids=ids, work_cache=str(tmp_path))
    assert got["acc"] == "GCF_1.1" and got["confidence"] == "validated"           # not the reference genome: the ids say so
    got = osdr.choose_genome(info, index=rows, ids=None, work_cache=str(tmp_path))
    assert got["acc"] == "GCF_2.1" and got["confidence"] in ("taxon", "species")   # no table: best metadata match
    got = osdr.choose_genome(info, index=rows, ids=["zzz1", "zzz2"], work_cache=str(tmp_path))
    assert got["confidence"] != "validated" and "no candidate matched" in got["evidence"]
    info["accessions"] = ["GCF_1.1"]
    assert osdr.choose_genome(info, index=rows, ids=None, work_cache=str(tmp_path))["confidence"] == "stated"


def test_choose_genome_prefers_better_named_candidate_within_margin(tmp_path, monkeypatch):
    rng = random.Random(8)
    seqs = [prot_seq(rng) for _ in range(10)]
    named = make_genome(tmp_path, "named", [(f"N{i}", "", f"g{i}", s) for i, s in enumerate(seqs[:8])])     # recognises 8/10
    other = make_genome(tmp_path, "other", [(f"N{i}", "", f"g{i}", s) for i, s in enumerate(seqs[:9])])     # recognises 9/10
    files = {"GCF_1.1": (tmp_path / "named.gff.gz", tmp_path / "named.faa.gz"), "GCF_2.1": (tmp_path / "other.gff.gz", tmp_path / "other.faa.gz")}
    monkeypatch.setattr(bacteria, "fetch_genome", lambda key, cache, g=None: files[key])
    rows = [row("GCF_1.1", "Foo bar", "S1", taxid="11", species_taxid="10"), row("GCF_2.1", "Foo bar", "", "reference genome", taxid="10", species_taxid="10")]
    info = dict(accessions=[], accession_digits=[], taxa={"10": "Foo bar"}, organisms=["Foo bar"], strains=["S1"])
    ids = [f"N{i}" for i in range(10)]
    got = osdr.choose_genome(info, index=rows, ids=ids, work_cache=str(tmp_path))
    assert got["acc"] == "GCF_1.1"                     # 80% vs 90%: within the margin, the name match wins (multi-strain arrays)
    got = osdr.choose_genome(info, index=rows, ids=ids[:2] + ["x"] * 8, work_cache=str(tmp_path))
    assert got["confidence"] != "validated"            # 20% matches nothing well enough
    only_other = ["N8"] * 6 + ["x"] * 2
    got = osdr.choose_genome(info, index=rows, ids=only_other, work_cache=str(tmp_path))
    assert got["acc"] == "GCF_2.1" and got["confidence"] == "validated"   # the named genome fails the threshold, the other doesn't


@pytest.mark.skipif(not shutil.which("diamond"), reason="diamond not installed")
def test_hub_only_build_cache_and_mapper(tmp_path, monkeypatch):
    rng = random.Random(11)
    base = [prot_seq(rng) for _ in range(6)]
    mk = lambda key, tag, drop: make_genome(tmp_path, key, [(f"{tag}{i}", "", f"g{i}", mutate(rng, s, 0.1)) for i, s in enumerate(base) if i != drop])  # noqa: E731
    reg = {k: dict(label=k.upper(), species=f"Sp {k}", acc=f"GCF_{n}.1", osds=[n]) for n, k in enumerate("abc", 1)}
    prots = {"a": mk("a", "A_", 5), "b": mk("b", "B_", 4), "c": mk("c", "C_", 3)}
    monkeypatch.setattr(bacteria, "fetch_genome", lambda key, cache, g=None: (tmp_path / f"{key}.gff.gz", tmp_path / f"{key}.faa.gz"))
    out, cache = tmp_path / "out", tmp_path / "cache"
    msgs = []
    bacteria.build_set(reg, str(out), str(cache), hub="a", hub_only=True, threads=1, progress=msgs.append, sensitivity="sensitive")
    assert sum("orthologs" in m for m in msgs) == 2                     # hub-only: n-1 comparisons
    assert not (out / "pairs" / "b__c.tsv").exists()
    m = bacteria.Mapper(str(out))
    assert m.hub == "a" and m.key_for("OSD-2") == "b" and m.key_for("C") == "c" and m.key_for(3) == "c"
    r = m.map(["B_0", "B_5", "nope"], "OSD-2", ["c"])                    # b -> c goes through the hub
    assert [x["locus_tag"] for x in r] == ["B_0", "B_5", ""]
    assert r[0]["hub_gene"] == "A_0" and r[0]["c"] == "C_0" and r[0]["c_pident"] > 50
    assert r[1]["hub_gene"] == "" and r[1]["c"] == ""                    # b's gene 5 has no hub ortholog
    assert (cache / "rbh").is_dir() and len(list((cache / "rbh").iterdir())) == 2
    calls = []
    real = bacteria.rbh
    monkeypatch.setattr(bacteria, "rbh", lambda *a, **k: calls.append(1) or real(*a, **k))
    bacteria.build_set(reg, str(tmp_path / "out2"), str(cache), hub="a", threads=1, sensitivity="sensitive")   # now all pairs
    assert len(calls) == 1                                                # only b-c is new; a-b and a-c come from the cache
    assert (tmp_path / "out2" / "pairs" / "b__c.tsv").exists()
