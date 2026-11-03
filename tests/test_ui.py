from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

PAGE = str(Path(__file__).resolve().parents[1] / "views" / "species_lookup.py")


def run_map(at, text):
    at.text_area[0].set_value(text).run()      # the Map button enables once there is input
    at.button[0].click().run(timeout=60)
    assert not at.exception, at.exception
    return at


def test_page_renders():
    at = AppTest.from_file(PAGE, default_timeout=60).run()
    assert not at.exception, at.exception
    assert at.title[0].value == "Look up genes across species"
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
        progress("[1/3] Loading F. bar S1 (GCF_000000001.1)")
        progress("    downloading F_genomic.gff.gz from https://ftp.ncbi.nlm.nih.gov/x/F_genomic.gff.gz")
        progress("[2/3] Comparing F. bar S1 with B. qux S2")
        progress("    reciprocal (each is the other's best hit): 3 ortholog pairs")
        progress("[3/3] Writing the result tables")
        hub_map = {"GCF_000000002.1": {h: (v[0], v[1]) for h, v in pair[("GCF_000000001.1", "GCF_000000002.1")].items()}}
        bacteria.write_tables(outdir, prot, hub_map, pair, "GCF_000000001.1", registry)

    monkeypatch.setattr(genomes, "assembly_index", lambda *a, **k: ROWS)
    monkeypatch.setattr(bacteria, "build_set", fake_build_set)
    return built


def test_any_genomes_page(stubs):
    at = AppTest.from_file(str(Path(PAGE).parent / "species_custom.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    at.text_area[0].set_value("GCF_000000001.1\nBaz qux S2\nnot a genome").run()
    assert [e.value for e in at.error] and "not a genome" in at.error[0].value
    assert [d["Genome"] for d in at.dataframe[0].value.to_dict("records")] == ["F. bar S1", "B. qux S2"]
    at.button[0].click().run()
    assert not at.exception, at.exception
    assert stubs == [(["GCF_000000001.1", "GCF_000000002.1"], "GCF_000000001.1")]
    assert any("in all 2 species" in s.value for s in at.subheader)             # results appear after the build
    log = [c.value for c in at.code if "Loading F. bar S1" in c.value]
    assert log and "downloading F_genomic.gff.gz from https://ftp.ncbi.nlm.nih.gov" in log[0] and "Finished in" in log[0]
    assert "[ " in log[0] or "s] " in log[0]                                    # every line is timestamped
    # the finished log must not vanish when the user keeps working (any later widget change reruns the page)
    at.text_area(key="any_text").set_value("F_0").run()
    assert not at.exception, at.exception
    assert any("Loading F. bar S1" in c.value and "Finished in" in c.value for c in at.code)


def test_any_osdr_page(stubs, monkeypatch):
    info = dict(organisms=["Foo bar"], strains=["S1"], tables=[("t.csv", "counts")])
    picks = {7: dict(acc="GCF_000000001.1", label="F. bar S1", species="Foo bar", confidence="validated", evidence="9/10 ids", info=info),
             8: dict(acc="GCF_000000002.1", label="B. qux S2", species="Baz qux", confidence="strain", evidence="name", info=dict(info, tables=[]))}
    def fake_lookup(studies, progress=None, **k):
        progress("[1/3] Loading the NCBI assembly index")
        progress("[2/3] OSD-7: reading the study and choosing its genome")
        progress("    fetching https://osdr.nasa.gov/osdr/data/osd/meta/7")
        return {k_: dict(v, osds=list(v["osds"])) for k_, v in REG.items()}, picks

    monkeypatch.setattr(osdr, "genomes_for", fake_lookup)
    at = AppTest.from_file(str(Path(PAGE).parent / "species_osdr.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    at.text_input[0].set_value("OSD-7, 8").run()
    at.button[0].click().run()
    assert not at.exception, at.exception
    rows = at.dataframe[0].value.to_dict("records")
    assert [r["Study"] for r in rows] == ["OSD-7", "OSD-8"] and rows[1]["Gene table"] == "raw reads only"
    lookup_log = lambda: [c.value for c in at.code if "fetching https://osdr.nasa.gov/osdr/data/osd/meta/7" in c.value]  # noqa: E731
    assert lookup_log() and "Finished in" in lookup_log()[0]
    at.checkbox(key="osdr_hubonly").check().run()                                # any later interaction reruns the page
    assert lookup_log(), "the look-up log must stay on the page"
    build = [b for b in at.button if b.label == "Build mapping"][0]
    build.click().run()
    assert not at.exception, at.exception
    assert stubs and any("in all 2 species" in s.value for s in at.subheader)
    at.text_area(key="osdr_any_text").set_value("F_0").run()
    assert any("Loading F. bar S1" in c.value for c in at.code)                 # build log still shown
    # and a changed study list hides the old logs (they belong to the old inputs)
    at.text_input[0].set_value("OSD-99").run()
    assert not any("Loading F. bar S1" in c.value for c in at.code)


# ---- explanation page (reads the committed tables in data/bacteria) ----

EXPL = str(Path(PAGE).parent / "species_explained.py")


def test_explanation_page_renders_with_live_examples():
    at = AppTest.from_file(EXPL, default_timeout=90).run()
    assert not at.exception, at.exception
    assert at.title[0].value == "How genes are matched between species"
    assert any("share the same symbol" in m.value for m in at.markdown)           # computed from the tables
    # default 'Try it' input is an OLD PA14 name: resolved to the current tag and mapped to B. subtilis
    assert any("is an old name for" in s.value for s in at.success)
    table = at.dataframe[-1].value
    assert table.loc[table["Genome"] == "B. subtilis 168", "Locus tag"].iloc[0] == "BSU_00010"


def test_explanation_shows_count_tables_and_protein_matching():
    at = AppTest.from_file(EXPL, default_timeout=90).run()
    assert not at.exception, at.exception
    frames = [d.value for d in at.dataframe]
    assert any(list(f.columns)[:2] == ["gene", "4.1"] and list(f["gene"])[0] == "PA14_RS00005" for f in frames)      # count table is a table
    assert any("trnS-Ser" in list(f.get("gene", [])) for f in frames)
    assert not any("PA14_RS00005\",186" in c.value for c in at.code)                                                # no comma-separated text
    text = " ".join(m.value for m in at.markdown)
    assert "Alignment: line the two proteins up" in text and "Reciprocal: ask the question back" in text
    assert any("picks BSU_00010, and BSU_00010 picks PA14_RS00005 back" in s.value for s in at.success)
    assert any("|" in c.value and "PA14" in c.value and "BSU" in c.value for c in at.code)                           # the real alignment


def test_explanation_try_it_unknown_and_probe():
    at = AppTest.from_file(EXPL, default_timeout=90).run()
    at.text_input[0].set_value("PA14_RS99999").run()
    assert not at.exception, at.exception
    assert any("not a gene in any" in w.value for w in at.warning)
    at.text_input[0].set_value("RruA0001").run()
    assert any("array probe, reporter or Entrez id" in s.value for s in at.success)


# ---- navigation and the between-strain pages ----

VIEWS = Path(PAGE).parent
ROOT = VIEWS.parent


def test_every_page_in_the_registry_exists_and_titles_are_unique():
    from gene_lingo import nav
    titles = [t for _, t, _ in nav.PAGES.values()]
    assert len(titles) == len(set(titles))
    for key, (path, title, icon) in nav.PAGES.items():
        assert (ROOT / path).exists(), path
    sectioned = [k for keys in nav.SECTIONS.values() for k in keys]
    assert sorted(sectioned) == sorted(nav.PAGES)                       # each page in exactly one section
    assert list(nav.SECTIONS) == ["Start", "Bacteria: between species", "Bacteria: between strains", "Human & animal genes"]


def test_router_starts_on_the_start_page():
    at = AppTest.from_file(str(ROOT / "app.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    assert at.title[0].value == "Gene Lingo"
    assert any("Two different questions" in h.value for h in at.header)


def test_home_numbers_come_from_the_tables():
    at = AppTest.from_file(str(VIEWS / "home.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    assert any("PAO1 and PA14 share 5,346 genes" in c.value for c in at.caption)


def test_strains_explained_page():
    at = AppTest.from_file(str(VIEWS / "strains_explained.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    assert at.title[0].value == "How genes are matched between strains"
    frames = [d.value for d in at.dataframe]
    summary = next(f for f in frames if "Only in A" in f.columns)
    assert len(summary) == 4 and summary.loc[0, "Shared"] == 5346 and summary.loc[0, "Only in B"] == 545
    assert any("T6SS" in " ".join(f.astype(str).values.ravel()) for f in frames)           # real PA14-only genes are shown


def test_compare_two_strains_page():
    at = AppTest.from_file(str(VIEWS / "strains_compare.py"), default_timeout=60).run()
    assert not at.exception, at.exception
    metrics = {m.label: m.value for m in at.metric}
    assert metrics["Shared"] == "5,346" and metrics["Only in A"] == "226" and metrics["Only in B"] == "545"
    at.text_area(key="strain_text").set_value("PA0001\nPA0002\nPA9999\nPA0050").run()
    at.button(key="strain_go").click().run()
    assert not at.exception, at.exception
    out = next(d.value for d in at.dataframe if "result" in d.value.columns)
    assert list(out["result"])[:3] == ["same gene in B", "same gene in B", "not recognised in A"]
    assert out.loc[0, "locus tag in B"] == "PA14_RS00005"
    at.selectbox[0].select("Escherichia coli").run()                                       # switching species keeps the page working
    assert not at.exception, at.exception
    assert {m.label: m.value for m in at.metric}["Shared"] == "3,889"
    at.selectbox[2].select("ecol_k12").run()                                                # same strain twice
    assert any("two different strains" in w.value for w in at.warning)
