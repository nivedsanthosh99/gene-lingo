"""Streamlit widgets shared by the bacterial pages: gene-mapping form, core/hub views, genome table."""
import csv
import io
import re

import pandas as pd
import streamlit as st

from . import bacteria


def read_names(upload, text):
    """Identifiers from an uploaded table (first column) or from pasted text."""
    if upload is not None:
        raw = upload.getvalue().decode("utf-8-sig")
        first = raw.splitlines()[0] if raw else ""
        delim = "\t" if "\t" in first else ";" if ";" in first else ","
        names = [r[0].strip().strip('"') for r in csv.reader(io.StringIO(raw), delimiter=delim) if r and r[0].strip()]
        if names and names[0].lower() in ("gene", "genes", "id", "locus_tag", "locus tag", "symbol", "gene_id", "input", "probeid", ""):
            names = names[1:]  # header row
        return names
    return [g for g in re.split(r"[\s,;]+", text) if g]


def label(mapper, k):
    g = mapper.genomes[k]
    studies = ", ".join(f"OSD-{o}" + ("*" if o in g.get("raw_only", []) else "") for o in g.get("osds", []))
    return f"{g['label']}  ({studies})" if studies else g["label"]


def hub_frame(mapper):
    """Hub table with one boolean column per species and `n_species` (including the hub's own)."""
    G, hub = mapper.genomes, mapper.hub
    df = pd.DataFrame(mapper.hub_table())
    others = [k for k in G if k != hub]
    for k in others:
        df[k] = df[k].fillna("")
    species = list(dict.fromkeys(G[k]["species"] for k in others))
    for sp in species:
        df[sp] = df[[k for k in others if G[k]["species"] == sp]].ne("").any(axis=1)
    df["n_species"] = df[species].sum(axis=1) + 1  # +1: the hub species itself
    return df, others, species


def results(mapper, key, about=""):
    """The four result tabs for a built set of tables."""
    G, hub = mapper.genomes, mapper.hub
    keys = list(G)
    hub_label = G[hub]["label"]
    hub_df, others, species = hub_frame(mapper)
    pretty = {k: G[k]["label"] for k in G}
    n_sp = len(species) + 1
    in_core = set(hub_df.loc[hub_df["n_species"] == n_sp, "locus_tag"])
    tab_map, tab_core, tab_hub, tab_about = st.tabs(["Map genes", "Core genes", f"Browse {hub_label} hub", "About"])

    with tab_map:
        c1, c2 = st.columns([1, 2])
        src_choice = c1.selectbox("Genes are from", ["Auto-detect", *keys], key=f"{key}_src",
                                  format_func=lambda k: k if k == "Auto-detect" else label(mapper, k))
        dsts = c2.multiselect("Map to", keys, default=[hub], key=f"{key}_dst", format_func=lambda k: label(mapper, k))
        if c2.checkbox("Map to all genomes", key=f"{key}_all"):
            dsts = keys
        text = st.text_area("Locus tags (one per line, or separated by commas/spaces)", height=140, key=f"{key}_text",
                            placeholder="One identifier per line: locus tags, old tags, array probe ids")
        upload = st.file_uploader("...or upload a CSV/TSV (e.g. a GeneLab counts table; first column is used)",
                                  type=["csv", "tsv", "txt"], key=f"{key}_up")
        names = read_names(upload, text)
        only_core = st.checkbox(f"Only genes in the {n_sp}-species core", key=f"{key}_core") if n_sp > 1 else False

        if st.button("Map", type="primary", disabled=not (names and dsts), key=f"{key}_go"):
            src = src_choice
            if src == "Auto-detect":
                src = mapper.detect(names)
                if src is None:
                    st.error("None of these identifiers match any of the genomes. Check the identifiers or pick a genome.")
                    st.stop()
                st.info(f"Detected **{G[src]['label']}**" + (f" ({label(mapper, src).split('(', 1)[1][:-1]})" if G[src].get("osds") else ""))
            df = pd.DataFrame(mapper.map(names, src, dsts))
            if only_core:
                df = df[df["hub_gene"].isin(in_core)]
            m1, m2, m3 = st.columns(3)
            m1.metric("Identifiers", len(set(names)))
            m2.metric("Recognised", int((df["locus_tag"] != "").sum()))
            m3.metric(f"Have a {hub_label} ortholog", int((df["hub_gene"] != "").sum()))
            st.write(" · ".join(f"**{G[d]['label']}** {int((df[d] != '').sum())}/{len(df)}" for d in dsts if d != src) + "  mapped")
            e = st.columns(3)
            show_id = e[0].checkbox("Show % identity", value=True, key=f"{key}_id")
            show_sym = e[1].checkbox("Show target symbols", value=False, key=f"{key}_sym")
            cols = ["input", "locus_tag", "symbol", "product", "hub_gene"]
            for d in dsts:
                cols += [d] + ([f"{d}_pident"] if show_id and d != src else []) + ([f"{d}_symbol"] if show_sym else [])
            names_ = {"hub_gene": f"{hub_label} hub gene"}
            for d in dsts:
                names_.update({d: G[d]["label"], f"{d}_pident": f"{G[d]['label']} %id", f"{d}_symbol": f"{G[d]['label']} symbol"})
            shown = df[list(dict.fromkeys(cols))].rename(columns=names_)
            tcols = [G[d]["label"] for d in dsts if G[d]["label"] in shown]
            st.dataframe(
                shown.style.format({c: "{:.1f}" for c in shown if c.endswith(" %id")}, na_rep="").apply(
                    lambda r: ["background-color: #fdd" if r["locus_tag"] == "" else
                               "background-color: #fff4d6" if c in tcols and r[c] == "" else "" for c in r.index], axis=1),
                width="stretch", hide_index=True)
            st.caption("Red: identifier not recognised. Yellow: recognised, no reciprocal best hit in that genome.")
            st.download_button("Download CSV", df.to_csv(index=False), "gene_lingo_bacteria.csv", "text/csv", key=f"{key}_dl")

    with tab_core:
        st.subheader(f"{len(in_core)} {hub_label} genes have an ortholog in all {n_sp} species")
        c1, c2 = st.columns(2)
        c1.markdown(f"**{hub_label} genes with an ortholog in each species**")
        c1.bar_chart(pd.Series({" ".join([sp.split()[0][0] + ".", sp.split()[1]]) if len(sp.split()) > 1 else sp: int(hub_df[sp].sum())
                                for sp in species}, name="genes"))
        c2.markdown(f"**Genes by number of species sharing them (incl. {hub_label.split()[0][0]}. {hub_label.split()[1]})**")
        c2.bar_chart(hub_df["n_species"].value_counts().sort_index().rename("genes"))
        q = st.text_input("Filter core genes (locus tag, symbol or product)", key=f"{key}_coreq")
        core_df = hub_df[hub_df["n_species"] == n_sp][["locus_tag", "symbol", "product", *others]]
        if q:
            core_df = core_df[core_df.apply(lambda r: q.lower() in " ".join(r.astype(str)).lower(), axis=1)]
        st.dataframe(core_df.rename(columns=pretty), width="stretch", hide_index=True)
        st.download_button("Download core genes (CSV)", core_df.to_csv(index=False), "core_genes.csv", "text/csv", key=f"{key}_coredl")

    with tab_hub:
        st.write(f"One row per {hub_label} gene; each column is the reciprocal-best-hit ortholog in that genome.")
        q = st.text_input("Search (locus tag, symbol or product)", key=f"{key}_hubq")
        view = hub_df[["locus_tag", "symbol", "product", *others, "n_species"]]
        if q:
            view = view[view.apply(lambda r: q.lower() in " ".join(r.astype(str)).lower(), axis=1)]
        st.caption(f"{len(view)} of {len(hub_df)} genes")
        st.dataframe(view.rename(columns=pretty), width="stretch", hide_index=True)
        st.download_button("Download table (CSV)", view.to_csv(index=False), "hub.csv", "text/csv", key=f"{key}_hubdl")

    with tab_about:
        if about:
            st.markdown(about)
        st.dataframe(pd.DataFrame([
            {"Genome": G[k]["label"], "Species": G[k]["species"], "Assembly": G[k].get("acc", ""),
             "OSD studies": ", ".join(f"OSD-{o}" for o in G[k].get("osds", [])), "Genes": len(mapper.genes[k].cds)} for k in G]),
            width="stretch", hide_index=True)


METHOD = """
**Method.** Protein sequences are compared with DIAMOND; two genes are orthologs when each is the other's best hit and
both alignments cover at least 50 % of both proteins (e-value ≤ 1e-5). Names are never compared. Old locus tags and
tags without underscores are resolved first. **Limits.** One ortholog per gene: paralog families can resolve to the
wrong member (check the % identity). Genes missing from the hub genome can still be mapped between other genomes but
are not in the hub table or the core.
"""


def build_and_show(registry, hub, hub_only, key, about="", cache=None):
    """'Build mapping' button for a genome registry; shows the results (and a zip of the tables) once built."""
    import json
    import os
    import tempfile
    import zipfile

    from . import genomes

    cache = cache or genomes.default_cache()
    sig = json.dumps([registry, hub, hub_only], sort_keys=True)
    n = len(registry)
    pairs = n - 1 if hub_only else n * (n - 1) // 2
    st.caption(f"{n} genomes, {pairs} ortholog comparisons (about {max(1, pairs * 3 // 60)}–{max(2, pairs * 6 // 60)} min the first time; "
               "genomes and finished comparisons are cached).")
    if st.button("Build mapping", type="primary", disabled=n < 2, key=f"{key}_build"):
        outdir = tempfile.mkdtemp(prefix="gene_lingo_")
        try:
            with st.status("Building ortholog tables...", expanded=True) as status:
                bacteria.build_set(registry, outdir, os.path.join(cache, "genomes"), hub, hub_only, 4, progress=status.write)
                status.update(label="Done", state="complete", expanded=False)
        except Exception as e:  # network / diamond / bad genome
            st.error(str(e))
            return
        st.session_state[f"{key}_built"] = (sig, outdir)
    built = st.session_state.get(f"{key}_built")
    if not (built and built[0] == sig):
        return
    outdir = built[1]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(outdir):
            for f in files:
                z.write(os.path.join(root, f), os.path.relpath(os.path.join(root, f), outdir))
    st.download_button("Download all tables (zip)", buf.getvalue(), "gene_lingo_tables.zip", "application/zip", key=f"{key}_zip")
    results(bacteria.Mapper(outdir), key, about)
