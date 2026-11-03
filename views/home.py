"""Start page: what the tool does, the two problems it solves, and which page to open."""
from pathlib import Path

import streamlit as st

from gene_lingo import bacteria, nav, strains

DATA = Path(__file__).resolve().parents[1] / "data" / "bacteria"

st.set_page_config(page_title="Gene Lingo", layout="centered")
st.title("Gene Lingo")
st.markdown(
    "Find the **same gene** in different bacteria, even when every dataset calls it something different. Built for "
    "NASA GeneLab (OSDR) bacterial studies; works for any bacterial genomes. Genes are matched by their **protein "
    "sequence**, never by name, so any locus-tag spelling or annotation works.")

M = bacteria.Mapper(str(DATA)) if (DATA / "genes.tsv.gz").exists() else None
pao1_pa14 = None
if M:
    c = strains.compare(M, "paer_pao1", "paer_pa14")
    pao1_pa14 = (len(c["shared"]), len(c["only_a"]), len(c["only_b"]), c["median_identity"])

st.header("Two different questions")
left, right = st.columns(2, gap="large")
with left:
    st.subheader("Between species")
    st.markdown("Is this gene of **B. subtilis** the same gene as one in **E. coli**?")
    st.markdown("Proteins are only about 40% identical, so the tool needs sensitive protein search and a common reference "
                "genome. Result: a table of the genes you can compare across all species.")
    st.caption("Example: *dnaA* is `BSU_00010`, `PA0001`, `PA14_RS00005`, `QV15_00005`, `b3702` and `RRU_RS00005`.")
with right:
    st.subheader("Between strains")
    st.markdown("One study used strain **PAO1**, another **PA14**. Which gene matches which, and which genes "
                "exist in only one?")
    st.markdown("Proteins are almost identical, so matching is easy; the work is knowing which strain a dataset used and "
                "spotting strain-specific genes.")
    if pao1_pa14:
        st.caption(f"Example: PAO1 and PA14 share {pao1_pa14[0]:,} genes (median {pao1_pa14[3]}% identical); "
                   f"{pao1_pa14[1]} are only in PAO1 and {pao1_pa14[2]} only in PA14.")

st.header("Which page do you need?")
rows = [
    ("I want to understand the problem and the method", "species_explained", "strains_explained"),
    ("I have gene IDs from an OSDR bacterial study (*B. subtilis*, *P. aeruginosa*, *S. aureus*, *E. coli*, *R. rubrum*)", "species_lookup", None),
    ("I have a different OSDR study and want its genome found automatically", "species_osdr", None),
    ("I have my own organisms or genomes", "species_custom", None),
    ("I want to compare two strains of one species", "strains_compare", None),
    ("My genes are human, mouse or another animal", "converter", None),
]
for text, key, key2 in rows:
    with st.container(border=True):
        st.markdown(text)
        links = st.columns(3)
        nav.link(links[0], key)
        if key2:
            nav.link(links[1], key2)

st.header("What is in the sidebar")
st.markdown(
    "- **Bacteria: between species** - how it works, then three tools (the five OSDR species, your own genomes, any OSDR study).\n"
    "- **Bacteria: between strains** - how it works, with a flowchart, then a tool to compare two strains.\n"
    "- **Human & animal genes** - the original gene-ID converter and ortholog finder.")
st.caption("Everything shown on the explanation pages comes from the project's real data tables.")
