from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

PAGE = str(Path(__file__).resolve().parents[1] / "pages" / "2_Bacterial_genomes.py")


def run_map(at, text):
    at.text_area[0].set_value(text).run()      # the Map button enables once there is input
    at.button[0].click().run(timeout=60)
    assert not at.exception, at.exception
    return at


def test_page_renders():
    at = AppTest.from_file(PAGE, default_timeout=60).run()
    assert not at.exception, at.exception
    assert at.title[0].value == "Bacterial genomes"
    assert any("in all 5 species" in s.value for s in at.subheader)


def test_map_old_and_new_names_autodetect():
    at = AppTest.from_file(PAGE, default_timeout=60).run()
    run_map(at, "PA14_00010\nPA14_RS00005\nbogus")
    assert any("Detected **P. aeruginosa PA14**" in i.value for i in at.info)
    df = at.dataframe[0].value
    assert list(df["locus_tag"]) == ["PA14_RS00005", "PA14_RS00005", ""]
    assert list(df["B. subtilis 168"]) == ["BSU_00010", "BSU_00010", ""]


def test_probe_ids_autodetect():
    at = AppTest.from_file(PAGE, default_timeout=60).run()
    run_map(at, "RruA0001")
    assert any("R. rubrum" in i.value for i in at.info)
    assert at.dataframe[0].value["locus_tag"][0] == "RRU_RS00005"


# ---- pages 3 and 4: any genomes / any OSDR study (network, assembly index and DIAMOND are stubbed) ----

import random  # noqa: E402

from helpers import make_genome, prot_seq  # noqa: E402

from gene_lingo import bacteria, genomes, osdr  # noqa: E402

ROWS = [dict(acc="GCF_000000001.1", category="reference genome", taxid="1", species_taxid="1", organism="Foo bar", strain="S1",
             level="Complete Genome", asm_name="x", date="2020", group="bacteria"),
        dict(acc="GCF_000000002.1", category="reference genome", taxid="2", species_taxid="2", organism="Baz qux", strain="S2",
             level="Complete Genome", asm_name="x", date="2020", group="bacteria")]
REG = {"GCF_000000001.1": dict(label="F. bar S1", species="Foo bar", acc="GCF_000000001.1", osds=[7]),
       "GCF_000000002.1": dict(label="B. qux S2", species="Baz qux", acc="GCF_000000002.1", osds=[8])}


@pytest.fixture
def stubs(tmp_path, monkeypatch):
    rng = random.Random(5)
    seqs = [prot_seq(rng) for _ in range(4)]
    prot = {"GCF_000000001.1": make_genome(tmp_path, "GCF_000000001.1", [(f"F_{i}", "", f"g{i}", s) for i, s in enumerate(seqs)]),
            "GCF_000000002.1": make_genome(tmp_path, "GCF_000000002.1", [(f"Q_{i}", "", f"g{i}", s) for i, s in enumerate(seqs[:3])])}
    pair = {("GCF_000000001.1", "GCF_000000002.1"): {f"F_{i}": (f"Q_{i}", 99.0) for i in range(3)}}
    built = []

    def fake_build_set(registry, outdir, cache, hub=None, hub_only=False, threads=4, progress=None, **kw):
        built.append((sorted(registry), hub))
        hub_map = {"GCF_000000002.1": {h: (v[0], v[1]) for h, v in pair[("GCF_000000001.1", "GCF_000000002.1")].items()}}
        bacteria.write_tables(outdir, prot, hub_map, pair, "GCF_000000001.1", registry)

    monkeypatch.setattr(genomes, "assembly_index", lambda *a, **k: ROWS)
    monkeypatch.setattr(bacteria, "build_set", fake_build_set)
    return built


def test_any_genomes_page(stubs):
    at = AppTest.from_file(str(Path(PAGE).parent / "3_Any_genomes.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    at.text_area[0].set_value("GCF_000000001.1\nBaz qux S2\nnot a genome").run()
    assert [e.value for e in at.error] and "not a genome" in at.error[0].value
    assert [d["Genome"] for d in at.dataframe[0].value.to_dict("records")] == ["F. bar S1", "B. qux S2"]
    at.button[0].click().run()
    assert not at.exception, at.exception
    assert stubs == [(["GCF_000000001.1", "GCF_000000002.1"], "GCF_000000001.1")]
    assert any("in all 2 species" in s.value for s in at.subheader)             # results appear after the build


def test_any_osdr_page(stubs, monkeypatch):
    info = dict(organisms=["Foo bar"], strains=["S1"], tables=[("t.csv", "counts")])
    picks = {7: dict(acc="GCF_000000001.1", label="F. bar S1", species="Foo bar", confidence="validated", evidence="9/10 ids", info=info),
             8: dict(acc="GCF_000000002.1", label="B. qux S2", species="Baz qux", confidence="strain", evidence="name", info=dict(info, tables=[]))}
    monkeypatch.setattr(osdr, "genomes_for", lambda studies, **k: ({k_: dict(v, osds=list(v["osds"])) for k_, v in REG.items()}, picks))
    at = AppTest.from_file(str(Path(PAGE).parent / "4_Any_OSDR_study.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    at.text_input[0].set_value("OSD-7, 8").run()
    at.button[0].click().run()
    assert not at.exception, at.exception
    rows = at.dataframe[0].value.to_dict("records")
    assert [r["Study"] for r in rows] == ["OSD-7", "OSD-8"] and rows[1]["Gene table"] == "raw reads only"
    build = [b for b in at.button if b.label == "Build mapping"][0]
    build.click().run()
    assert not at.exception, at.exception
    assert stubs and any("in all 2 species" in s.value for s in at.subheader)
