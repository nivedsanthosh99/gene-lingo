"""Streamlit page: one-page summary of the project: problem, method, results, status and open work."""
import csv
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import streamlit as st

from gene_lingo import bacteria, nav, strains
from gene_lingo.ui import hub_frame

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "bacteria"
ALL = ROOT / "data" / "osdr_all"
SURVEY = ROOT / "data" / "osdr_bacterial_survey.tsv"

st.set_page_config(page_title="Project summary - Gene Lingo", layout="centered")
st.title("Project summary")
st.markdown(
    "**Gene Lingo** matches the same gene across bacterial species and strains, so that gene-activity results from "
    "different NASA GeneLab (OSDR) experiments can be compared. Datasets name genes differently (different species, "
    "strains and annotation versions), so the tool matches genes by **protein sequence** and reports every result in the "
    "locus tags the datasets actually use.")

# ---- live numbers
have = (DATA / "genes.tsv.gz").exists()
if have:
    M = bacteria.Mapper(str(DATA))
    hub_df, others, species = hub_frame(M)
    core_n = int((hub_df["n_species"] == len(species) + 1).sum())
    pairs = strains.all_pairs(M)
srv = []
if SURVEY.exists():
    with open(SURVEY) as fh:
        srv = list(csv.DictReader(fh, delimiter="\t"))
conf = Counter(r["confidence"] for r in srv)
n_all = len(json.loads((ALL / "genomes.json").read_text())["genomes"]) if (ALL / "genomes.json").exists() else None

st.header("At a glance")
t = st.columns(4)
if have:
    t[0].metric("Species / genomes", f"{len(species) + 1} / {len(M.genomes)}")
    t[1].metric("Core genes (all 5 species)", f"{core_n}")
t[2].metric("OSDR studies surveyed", f"{len(srv)}")
if n_all:
    t[3].metric("Genomes across them", f"{n_all}")
st.caption("Core genes are *B. subtilis* genes that have a counterpart in every one of the five species. "
           "Studies are all bacterial transcription studies in OSDR.")

# ---- problem / approach
st.header("The problem")
st.markdown(
    "- The same gene has a different name in every species: *dnaA* is `BSU_00010`, `PA0001`, `PA14_RS00005`, `QV15_00005`, "
    "`b3702`, `RRU_RS00005`.\n"
    "- One genome can have several spellings (`PA14_RS00005` and the older `PA14_00010`), so tools keyed to one annotation "
    "(such as Ortholuge) miss genes.\n"
    "- Gene symbols disagree between species (`metS` / `metG`), so names cannot be the key.\n"
    "- Studies use different strains, some have no gene table, and arrays report probe IDs instead of locus tags.")
st.header("The approach")
st.markdown(
    "1. **Find each dataset's genome** from the OSDR metadata, and check it against the identifiers in the study's table.\n"
    "2. **Build a protein list** per genome in which every protein carries the dataset's own locus tag (old names and probe IDs resolve).\n"
    "3. **Match proteins** between genomes by reciprocal best hit (DIAMOND, e-value 1e-5, at least 50% of both proteins covered).\n"
    "4. **Use *B. subtilis* as a hub** so all species can be written in one set of tags and a core gene set defined; pairs are also compared directly.\n"
    "5. For **strains of one species** the same matching shows which genes are shared and which exist in only one strain.")

st.header("Inputs, tools and outputs")
st.dataframe(pd.DataFrame([
    {"Inputs": "OSD study numbers, or genome accessions / organism names / your own annotation files; gene IDs to translate",
     "Tools": "Python, DIAMOND (protein search), NCBI genomes and RefSeq assembly index, OSDR API, Streamlit",
     "Outputs": "hub table, core genes, a table per genome pair, ID translator, strain comparisons, a survey of OSDR studies"}]),
    hide_index=True, width="stretch")
st.markdown(
    "**General or OSDR-specific?** Both. The matching engine and genome lookup work on any bacterial or archaeal genomes. "
    "The OSDR layer (study metadata, genome inference, the array-probe rules for OSD-14, 15, 31, 39 and 43) is specific to OSDR.")

# ---- results
st.header("What it found")
if have:
    hub_label = M.genomes[M.hub]["label"]
    counts = {M.genomes[k]["label"]: int(hub_df[k].ne("").sum()) for k in others}
    st.markdown(f"**Between species.** Each other genome has an ortholog for these many of the {len(M.genes[M.hub].cds):,} {hub_label} genes:")
    st.dataframe(pd.DataFrame({"Genome": list(counts), "Genes with an ortholog": list(counts.values())}), hide_index=True, width="stretch")
    pa = next(r for r in pairs if r["strain_a"].endswith("PAO1"))
    st.markdown(
        f"**Between strains.** *P. aeruginosa* PAO1 and PA14 share {pa['shared']:,} genes (median {pa['median_identity']:.1f}% identical); "
        f"{pa['only_a']} are only in PAO1 and {pa['only_b']} only in PA14. The three *E. coli* strains share 3,766 to 3,889 genes "
        "per pair at a median of 100% identity.")
if srv:
    st.markdown(
        f"**OSDR survey.** Genome chosen for {len(srv)} bacterial studies: {conf.get('validated', 0)} validated against the study's own "
        f"identifiers, {conf.get('stated', 0)} named in the study description, {conf.get('strain', 0)} by strain name and "
        f"{conf.get('species', 0) + conf.get('taxon', 0)} at species level. The unvalidated ones have no gene table to check.")

# ---- validation
st.header("Checks done, and what is not checked yet")
c1, c2 = st.columns(2)
with c1:
    st.markdown("**Done**")
    st.markdown(
        "- 92-98% of identifiers in each study's table are genes of the chosen genome (the rest are tRNA/rRNA/pseudogenes)\n"
        "- Strains of one species match at 99.6-100% median identity\n"
        "- Through the hub or directly gives the same gene for 94.7% of genes\n"
        "- 96.6% of the 3,786 *R. rubrum* array oligos lie inside the gene assigned by the probe rule\n"
        "- Well-known conserved genes (*rpoB*, *recA*, *gyrA*, *groEL*) are in the core\n"
        "- Automated tests for the engine, the OSDR logic and every page")
with c2:
    st.markdown("**Not done yet**")
    st.markdown(
        "- No comparison with Ortholuge or OrthoFinder, and the cut-offs are conventional, not tuned\n"
        "- No test of the strain/genome inference with the metadata hidden\n"
        "- No check of gene neighbourhoods (true orthologs usually sit next to each other)\n"
        "- The genome for studies without a gene table is a name-based guess\n"
        "- Not yet applied to real expression tables")

st.header("Limits")
st.markdown(
    "- One best match per gene, so gene families can resolve to the wrong member (*lysS* matched to *lysU*).\n"
    "- Very diverged genes can be missed (*R. rubrum* *ftsZ* is not in the hub table).\n"
    "- Only protein-coding genes are mapped.\n"
    "- \"Only in one strain\" means no counterpart was found, not that the gene is absent.")

st.header("Status and next steps")
st.dataframe(pd.DataFrame([
    {"Area": "Matching engine, tables, OSDR genome inference", "Status": "Working and tested"},
    {"Area": "Web app (explanations, lookup, your own genomes, OSDR studies, strains)", "Status": "Working and tested; hosted on a free Streamlit or Hugging Face host"},
    {"Area": "Validation against other tools", "Status": "Not done: needs your Ortholuge / existing mappings, or OrthoFinder"},
    {"Area": "Translate a whole counts / differential-expression table and merge studies", "Status": "Not built (the main missing step for research use)"},
    {"Area": "Other OSDR data (genome sequencing, proteomics, plants, mouse, human)", "Status": "Not built"}]),
    hide_index=True, width="stretch")

st.header("Go to")
cols = st.columns(3)
nav.link(cols[0], "species_explained")
nav.link(cols[1], "strains_explained")
nav.link(cols[2], "species_osdr")
st.caption("More detail: `docs/METHOD.md` (technical method) and `README.md` (how to run and deploy) in the repository.")
