"""Streamlit UI: streamlit run app.py"""
import csv
import io
import re

import pandas as pd
import streamlit as st

from gene_lingo import ID_TYPES
from gene_lingo.pipeline import translate
from gene_lingo.species import SPECIES

st.set_page_config(page_title="gene-lingo", layout="wide")
st.title("gene-lingo")
st.caption("Convert gene IDs and find orthologs. Same input and output species = ID conversion only.")

names = [s[0] for s in SPECIES]
c1, c2, c3, c4 = st.columns(4)
in_species = c1.selectbox("Input species", names, index=names.index("human"))
out_species = c2.selectbox("Output species", names, index=names.index("human"))
in_type = c3.selectbox("Input ID type", ["auto", *ID_TYPES], index=0)
out_type = c4.selectbox("Output ID type", list(ID_TYPES), index=0)
ortho_type = "all"
if in_species != out_species:
    ortho_type = st.selectbox("Orthology type", ["all", "one2one", "one2many", "many2many"])

text = st.text_area("Genes (one per line, or separated by commas/spaces)", height=150)
upload = st.file_uploader("...or upload a CSV/text file (first column is used)", type=["csv", "txt", "tsv"])

genes = []
if upload is not None:
    raw = upload.getvalue().decode("utf-8-sig")
    delim = "\t" if "\t" in raw.splitlines()[0] else ";" if ";" in raw.splitlines()[0] else ","
    reader = csv.reader(io.StringIO(raw), delimiter=delim)
    genes = [r[0].strip() for r in reader if r and r[0].strip()]
    if genes and genes[0].lower() in ("gene", "genes", "id", "symbol", "gene_id", "input"):
        genes = genes[1:]  # header row
else:
    genes = [g for g in re.split(r"[\s,;]+", text) if g]

if st.button("Run", type="primary", disabled=not genes):
    try:
        with st.spinner(f"Translating {len(genes)} genes..."):
            df = pd.DataFrame(translate(genes, in_species, out_species, in_type, out_type, ortho_type))
    except Exception as e:  # network / API / bad combination
        st.error(str(e))
    else:
        mapped = df.loc[df["output_id"] != "", "input"].nunique()
        st.write(f"**{mapped} of {df['input'].nunique()}** genes mapped")
        st.dataframe(df.style.apply(
            lambda r: ["background-color: #fdd" if r["output_id"] == "" else "" for _ in r], axis=1),
            use_container_width=True, hide_index=True)
        st.download_button("Download CSV", df.to_csv(index=False), "gene_lingo_results.csv", "text/csv")
