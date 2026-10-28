"""Streamlit page: map OSDR bacterial locus tags between genomes (tables built by `gene-lingo bacteria build`)."""
from pathlib import Path

import streamlit as st

from gene_lingo import bacteria
from gene_lingo.ui import METHOD, results

DATA = Path(__file__).resolve().parents[1] / "data" / "bacteria"

st.set_page_config(page_title="Bacterial genomes - gene-lingo", layout="wide")
st.title("Bacterial genomes")
st.caption("Map locus tags between the genomes behind the OSDR bacterial studies. Orthologs are reciprocal best "
           "protein hits, so old names, other annotations and array probe ids all work. For other genomes see "
           "*Any genomes*; for other OSDR studies see *Any OSDR study*.")

if not (DATA / "genes.tsv.gz").exists():
    st.error(f"Mapping tables not found in `{DATA}`. Build them with `gene-lingo bacteria build -o data/bacteria`.")
    st.stop()


@st.cache_resource
def get_mapper():
    return bacteria.Mapper(str(DATA))


results(get_mapper(), "osdr8", about=METHOD + """
Array probe ids (`PA0001_dnaA_at`, `RruA0001`) and Entrez ids (OSD-43) are resolved too. Studies marked * have no
processed gene table on OSDR; their genome is inferred from the study metadata.
""")
