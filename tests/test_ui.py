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
