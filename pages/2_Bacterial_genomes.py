"""Streamlit page: map OSDR bacterial locus tags between genomes (tables built by `gene-lingo bacteria build`)."""
import csv
import io
import re
from pathlib import Path

import pandas as pd
import streamlit as st

from gene_lingo import bacteria

DATA = Path(__file__).resolve().parents[1] / "data" / "bacteria"
RAW_ONLY = {215, 216, 728}  # no processed gene tables on OSDR; genome inferred from the study metadata

st.set_page_config(page_title="Bacterial genomes - gene-lingo", layout="wide")
st.title("Bacterial genomes")
st.caption("Map locus tags between the genomes behind the OSDR bacterial studies. Orthologs are reciprocal best "
           "protein hits, so old names, other annotations and array probe ids all work.")

if not (DATA / "genes.tsv.gz").exists():
    st.error(f"Mapping tables not found in `{DATA}`. Build them with `gene-lingo bacteria build -o data/bacteria`.")
    st.stop()


@st.cache_resource
def get_mapper():
    return bacteria.Mapper(str(DATA))


mapper = get_mapper()
G = bacteria.GENOMES
HUB = bacteria.HUB


def osd_text(k):
    return ", ".join(f"OSD-{o}" + ("*" if o in RAW_ONLY else "") for o in G[k]["osds"])


def label(k):
    return f"{G[k]['label']}  ({osd_text(k)})"


def read_names(upload, text):
    if upload is not None:
        raw = upload.getvalue().decode("utf-8-sig")
        first = raw.splitlines()[0] if raw else ""
        delim = "\t" if "\t" in first else ";" if ";" in first else ","
        names = [r[0].strip().strip('"') for r in csv.reader(io.StringIO(raw), delimiter=delim) if r and r[0].strip()]
        if names and names[0].lower() in ("gene", "genes", "id", "locus_tag", "symbol", "gene_id", "input", "probeid", ""):
            names = names[1:]  # header row
        return names
    return [g for g in re.split(r"[\s,;]+", text) if g]


hub_rows = mapper.hub_table()
hub_df = pd.DataFrame(hub_rows).rename(columns={f"{HUB}_locus_tag": "locus_tag"})
others = [k for k in G if k != HUB]
species_of = {k: G[k]["species"] for k in others}
for k in others:
    hub_df[k] = hub_df[k].fillna("")
species = list(dict.fromkeys(species_of.values()))
pretty = {k: G[k]["label"] for k in G}
for sp in species:
    hub_df[sp] = hub_df[[k for k in others if species_of[k] == sp]].ne("").any(axis=1)
hub_df["n_species"] = hub_df[species].sum(axis=1) + 1  # +1: B. subtilis itself
in_core = set(hub_df.loc[hub_df["n_species"] == len(species) + 1, "locus_tag"])

tab_map, tab_core, tab_hub, tab_about = st.tabs(["Map genes", "Core genes", "Browse B. subtilis hub", "About"])

# ---------------------------------------------------------------- map
with tab_map:
    keys = list(G)
    c1, c2 = st.columns([1, 2])
    src_choice = c1.selectbox("Genes are from", ["Auto-detect", *keys], format_func=lambda k: k if k == "Auto-detect" else label(k))
    dsts = c2.multiselect("Map to", keys, default=[HUB], format_func=label)
    if c2.checkbox("Map to all genomes"):
        dsts = keys
    text = st.text_area("Locus tags (one per line, or separated by commas/spaces)", height=140,
                        placeholder="PA14_RS00005\nPA14_00010\nPA0001_dnaA_at\nBSU00010\nRruA0001")
    upload = st.file_uploader("...or upload a CSV/TSV (e.g. a GeneLab counts table; first column is used)",
                              type=["csv", "tsv", "txt"])
    names = read_names(upload, text)
    only_core = st.checkbox("Only genes in the 5-species core", value=False)

    if st.button("Map", type="primary", disabled=not (names and dsts)):
        src = src_choice
        if src == "Auto-detect":
            src = bacteria.detect_genome(names, mapper.genes)
            if src is None:
                st.error("None of these identifiers match any of the genomes. Check the identifiers or pick a genome.")
                st.stop()
            st.info(f"Detected **{G[src]['label']}** ({osd_text(src)})")
        rows = mapper.map(names, src, dsts)
        df = pd.DataFrame(rows)
        if only_core:
            df = df[df["bsub_168_hub"].isin(in_core)]
        resolved = int((df["locus_tag"] != "").sum())
        m1, m2, m3 = st.columns(3)
        m1.metric("Identifiers", len(set(names)))
        m2.metric("Recognised", f"{resolved}")
        m3.metric("Have a B. subtilis ortholog", int((df["bsub_168_hub"] != "").sum()))
        st.write(" · ".join(f"**{G[d]['label']}** {int((df[d] != '').sum())}/{len(df)}" for d in dsts if d != src)
                 + "  mapped")
        extras = st.columns(3)
        show_id = extras[0].checkbox("Show % identity", value=True)
        show_sym = extras[1].checkbox("Show target symbols", value=False)
        cols = ["input", "locus_tag", "symbol", "product", "bsub_168_hub"]
        for d in dsts:
            cols += [d] + ([f"{d}_pident"] if show_id and d != src else []) + ([f"{d}_symbol"] if show_sym else [])
        names_ = {"bsub_168_hub": "B. subtilis hub gene"}
        for d in dsts:
            names_.update({d: G[d]["label"], f"{d}_pident": f"{G[d]['label']} %id", f"{d}_symbol": f"{G[d]['label']} symbol"})
        shown = df[list(dict.fromkeys(cols))].rename(columns=names_)
        tcols = [G[d]["label"] for d in dsts if G[d]["label"] in shown]
        st.dataframe(
            shown.style.format({c: "{:.1f}" for c in shown if c.endswith(" %id")}, na_rep="").apply(lambda r: ["background-color: #fdd" if r["locus_tag"] == "" else
                                         "background-color: #fff4d6" if c in tcols and r[c] == "" else ""
                                         for c in r.index], axis=1),
            width="stretch", hide_index=True)
        st.caption("Red: identifier not recognised. Yellow: recognised, no reciprocal best hit in that genome.")
        st.download_button("Download CSV", df.to_csv(index=False), "gene_lingo_bacteria.csv", "text/csv")

# ---------------------------------------------------------------- core
with tab_core:
    st.subheader(f"{len(in_core)} B. subtilis genes have an ortholog in all {len(species) + 1} species")
    c1, c2 = st.columns(2)
    c1.markdown("**B. subtilis genes with an ortholog in each species**")
    c1.bar_chart(pd.Series({" ".join([sp.split()[0][0] + ".", sp.split()[1]]): int(hub_df[sp].sum()) for sp in species}, name="genes"))
    c2.markdown("**Genes by number of species sharing them (incl. B. subtilis)**")
    c2.bar_chart(hub_df["n_species"].value_counts().sort_index().rename("genes"))
    q = st.text_input("Filter core genes (locus tag, symbol or product)", key="core_q")
    core_df = hub_df[hub_df["n_species"] == len(species) + 1][["locus_tag", "symbol", "product", *others]]
    if q:
        core_df = core_df[core_df.apply(lambda r: q.lower() in " ".join(r.astype(str)).lower(), axis=1)]
    st.dataframe(core_df.rename(columns=pretty), width="stretch", hide_index=True)
    st.download_button("Download core genes (CSV)", core_df.to_csv(index=False), "core_genes.csv", "text/csv")

# ---------------------------------------------------------------- hub
with tab_hub:
    st.write("One row per *B. subtilis* 168 gene; each column is the reciprocal-best-hit ortholog in that genome.")
    q = st.text_input("Search (locus tag, symbol or product)", key="hub_q")
    view = hub_df[["locus_tag", "symbol", "product", *others, "n_species"]]
    if q:
        view = view[view.apply(lambda r: q.lower() in " ".join(r.astype(str)).lower(), axis=1)]
    st.caption(f"{len(view)} of {len(hub_df)} genes")
    st.dataframe(view.rename(columns=pretty), width="stretch", hide_index=True)
    st.download_button("Download table (CSV)", view.to_csv(index=False), "bsub_hub.csv", "text/csv")

# ---------------------------------------------------------------- about
with tab_about:
    st.markdown("""
**Method.** Each OSD is tied to the genome whose locus tags appear in its processed data (NCBI GFF + protein FASTA).
Protein sequences are compared with DIAMOND; two genes are orthologs when each is the other's best hit and both
alignments cover at least 50 % of both proteins (e-value ≤ 1e-5). This is done for every pair of genomes and against
the hub genome *B. subtilis* 168. Names are never compared. Old locus tags (`PA14_00010`), tags without underscores
(`BSU00010`), array probe ids (`PA0001_dnaA_at`, `RruA0001`) and Entrez ids (OSD-43) are resolved first.

**Limits.** One ortholog per gene: paralog families can resolve to the wrong member (check the % identity).
Genes missing from *B. subtilis* can still be mapped between other genomes (use the pair tables), but are not in the hub
or the core. Studies marked * have no processed gene table on OSDR; their genome is inferred from the study metadata.
""")
    st.dataframe(pd.DataFrame([
        {"Genome": G[k]["label"], "Species": G[k]["species"], "Assembly": G[k]["acc"], "OSD studies": osd_text(k),
         "Genes": len(mapper.genes[k].cds)} for k in G]), width="stretch", hide_index=True)
