"""Streamlit page: compare two strains of one species: shared genes, strain-specific genes, ID translator."""
from pathlib import Path

import pandas as pd
import streamlit as st

from gene_lingo import bacteria, nav, strains
from gene_lingo.ui import read_names

DATA = Path(__file__).resolve().parents[1] / "data" / "bacteria"

st.set_page_config(page_title="Compare two strains - Gene Lingo", layout="wide")
st.title("Compare two strains")
st.caption("Pick two strains of one species. See which genes they share, which exist in only one, and translate gene IDs "
           "from one strain to the other. How it works: ")
nav.link(st, "strains_explained")

if not (DATA / "genes.tsv.gz").exists():
    st.error(f"Mapping tables not found in `{DATA}`. Build them with `gene-lingo bacteria build -o data/bacteria`.")
    st.stop()


@st.cache_resource
def get_mapper():
    return bacteria.Mapper(str(DATA))


M = get_mapper()
G = M.genomes
groups = strains.species_groups(M)
if not groups:
    st.warning("The loaded tables contain only one strain per species.")
    st.stop()


def label(k):
    studies = ", ".join(f"OSD-{o}" for o in G[k].get("osds", []))
    return f"{G[k]['label']}  ({studies})" if studies else G[k]["label"]


c1, c2, c3 = st.columns(3)
species = c1.selectbox("Species", list(groups))
keys = groups[species]
a = c2.selectbox("Strain A", keys, index=0, format_func=label)
b = c3.selectbox("Strain B", keys, index=1, format_func=label)
if a == b:
    st.warning("Choose two different strains.")
    st.stop()

cmp = strains.compare(M, a, b)
A, B = G[a]["label"], G[b]["label"]
m = st.columns(6)
m[0].metric(f"Genes in {A.split()[-1]}", f"{cmp['n_a']:,}")
m[1].metric(f"Genes in {B.split()[-1]}", f"{cmp['n_b']:,}")
m[2].metric("Shared", f"{len(cmp['shared']):,}")
m[3].metric("Only in A", f"{len(cmp['only_a']):,}")
m[4].metric("Only in B", f"{len(cmp['only_b']):,}")
m[5].metric("Median % identical", f"{cmp['median_identity']:.1f}")
st.caption(f"A = {A}, B = {B}. \"Only in\" means no counterpart was found by reciprocal best hit, which is not proof that the "
           "gene is absent (it may be too diverged or a pseudogene).")

left, right = st.columns([2, 3])
with left:
    st.markdown("**How identical are the shared genes?**")
    hist = pd.DataFrame(strains.identity_histogram(cmp["shared"]), columns=["% identical", "genes"])
    st.bar_chart(hist, x="% identical", y="genes", sort=False, horizontal=True)

tab_tr, tab_sh, tab_a, tab_b = st.tabs(["Translate IDs", "Shared genes", f"Only in A", f"Only in B"])


with tab_tr:
    text = st.text_area(f"Locus tags of strain A ({A}), one per line, or upload a table below. Old names and probe IDs work too.",
                        height=130, key="strain_text")
    upload = st.file_uploader("...or upload a CSV/TSV (first column is used)", type=["csv", "tsv", "txt"], key="strain_up")
    names = read_names(upload, text)
    if st.button("Translate A → B", type="primary", disabled=not names, key="strain_go"):
        rows = M.map(names, a, [b])
        df = pd.DataFrame(rows)
        status = []
        for r in rows:
            status.append("not recognised in A" if not r["locus_tag"] else f"same gene in B" if r[b] else "only in A (no counterpart in B)")
        out = pd.DataFrame({"input": df["input"], "locus tag in A": df["locus_tag"], "symbol": df["symbol"], "product": df["product"],
                            "locus tag in B": df[b], "% identical": df[f"{b}_pident"], "result": status})
        st.write(f"**{sum(s == 'same gene in B' for s in status)}** of {len(out)} map to B; "
                 f"{sum(s.startswith('only') for s in status)} only in A; {sum(s.startswith('not') for s in status)} not recognised.")
        st.dataframe(out.style.format({"% identical": "{:.1f}"}, na_rep="").apply(
            lambda r: ["background-color: #fdd" if r["result"].startswith("not") else
                       "background-color: #fff4d6" if r["result"].startswith("only") else "" for _ in r], axis=1),
            hide_index=True, width="stretch")
        st.download_button("Download CSV", out.to_csv(index=False), "strain_translation.csv", "text/csv", key="strain_dl")

with tab_sh:
    sh = pd.DataFrame([{"locus tag in A": ta, "locus tag in B": tb, "symbol": M.genes[a].cds[ta]["symbol"],
                        "product": M.genes[a].cds[ta]["product"], "% identical": pid} for ta, tb, pid in cmp["shared"]])
    q = st.text_input("Filter (locus tag, symbol or product)", key="strain_shq")
    if q:
        sh = sh[sh.apply(lambda r: q.lower() in " ".join(r.astype(str)).lower(), axis=1)]
    st.caption(f"{len(sh):,} of {len(cmp['shared']):,} shared genes")
    st.dataframe(sh.style.format({"% identical": "{:.1f}"}), hide_index=True, width="stretch")
    st.download_button("Download shared genes (CSV)", sh.to_csv(index=False), "shared_genes.csv", "text/csv", key="strain_shdl")

for tab, tags, key, name, kk in ((tab_a, cmp["only_a"], a, A, "a"), (tab_b, cmp["only_b"], b, B, "b")):
    with tab:
        df = pd.DataFrame([{"locus tag": t, "symbol": M.genes[key].cds[t]["symbol"], "product": M.genes[key].cds[t]["product"]} for t in tags])
        hide = st.checkbox('Hide "hypothetical protein"', value=False, key=f"strain_hyp_{kk}")
        if hide and not df.empty:
            df = df[~df["product"].str.contains("hypothetical", case=False)]
        q = st.text_input("Filter (locus tag, symbol or product)", key=f"strain_only_q_{kk}")
        if q and not df.empty:
            df = df[df.apply(lambda r: q.lower() in " ".join(r.astype(str)).lower(), axis=1)]
        st.caption(f"{len(df):,} of {len(tags):,} genes found only in {name}")
        st.dataframe(df, hide_index=True, width="stretch")
        st.download_button(f"Download (CSV)", df.to_csv(index=False), f"only_in_{kk}.csv", "text/csv", key=f"strain_only_dl_{kk}")
