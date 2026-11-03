"""Streamlit page: the between-strain problem, a flowchart of the method, and live examples."""
from pathlib import Path

import pandas as pd
import streamlit as st

from gene_lingo import bacteria, nav, strains

DATA = Path(__file__).resolve().parents[1] / "data" / "bacteria"

st.set_page_config(page_title="How strain matching works - Gene Lingo", layout="centered")

if not (DATA / "genes.tsv.gz").exists():
    st.title("How genes are matched between strains")
    st.error(f"Mapping tables not found in `{DATA}`. Build them with `gene-lingo bacteria build -o data/bacteria`.")
    st.stop()


@st.cache_resource
def get_mapper():
    return bacteria.Mapper(str(DATA))


M = get_mapper()
G = M.genomes
groups = strains.species_groups(M)
COLOUR = {"Pseudomonas aeruginosa": "blue", "Escherichia coli": "violet"}


def badge(tag, species):
    return f":{COLOUR.get(species, 'gray')}-badge[{tag.replace('_', chr(92) + '_')}]"


st.title("How genes are matched between strains")
st.markdown(
    "A **strain** is one variety of a species, for example *P. aeruginosa* PAO1 and PA14, or *E. coli* K-12, B and ATCC 4157. "
    "Different GeneLab studies often use different strains of the same species, and each strain's genome is annotated "
    "with its own locus tags. To combine their results you need to know which gene in one strain is the same gene in the "
    "other.")

# ---- the problem
st.header("The problem")
st.write("Which strain each study used (from the OSDR metadata):")
st.markdown(
    "| Species | Strain | OSDR studies | Locus tags look like |\n|---|---|---|---|\n"
    f"| *P. aeruginosa* | PAO1 | OSD-14, OSD-15 | {badge('PA0001', 'Pseudomonas aeruginosa')} |\n"
    f"| *P. aeruginosa* | PA14 | OSD-554 | {badge('PA14_RS00005', 'Pseudomonas aeruginosa')} |\n"
    f"| *E. coli* | K-12 MG1655 | OSD-43 (and OSD-452 by name) | {badge('b3702', 'Escherichia coli')} |\n"
    f"| *E. coli* | ATCC 4157 (NCTC 86) | OSD-95 | {badge('E4V11_RS00005', 'Escherichia coli')} |\n"
    f"| *E. coli* | B REL606 | OSD-728 | {badge('ECB_RS18960', 'Escherichia coli')} |")
st.write("Those are all the same gene, *dnaA*, written once per strain. The same strain can also have several spellings "
         "(PA14's *dnaA* is `PA14_RS00005` today and was `PA14_00010` before), and some genes exist in one strain only.")

st.subheader("How this differs from the between-species problem")
st.markdown(
    "| | Between species | Between strains |\n|---|---|---|\n"
    "| How similar are the proteins? | median about 40% identical (most between 27% and 70%) | usually above 99% identical |\n"
    "| Hard part | finding the counterpart among distant proteins | knowing which strain a study used, and genes found in only one strain |\n"
    "| Typical result | about 950-1,400 of 4,237 *B. subtilis* genes have a counterpart | most genes match (PAO1 and PA14: 96% of PAO1's genes) |\n"
    "| Common reference genome needed? | yes (*B. subtilis*, to define a core) | no, strains are compared directly |")

# ---- flowchart
st.header("The method")
st.graphviz_chart(r'''
digraph G {
  rankdir=TB; bgcolor="transparent"; nodesep=0.3; ranksep=0.38;
  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=11, fillcolor="#eef2f7", color="#6b7c93", fontcolor="#15262c"];
  edge [fontname="Helvetica", fontsize=10, color="#6b7c93"];
  s  [label="Study with strain A's identifiers\ne.g. OSD-554: PA14_RS00005, PA14_00010", fillcolor="#dff0e3"];
  w  [label="Which strain was it really?\nstrain named in the OSDR metadata,\nthen checked: do the study's IDs exist in that genome?"];
  ok [label="Enough IDs found?\n(table available, >= 45% recognised)", shape=diamond, fillcolor="#fff4d6"];
  v  [label="Strain confirmed"];
  g  [label="Strain taken from the metadata only\n(raw-reads studies cannot be checked)", fillcolor="#fff4d6"];
  r  [label="Resolve spelling\nold names, current names, probe IDs"];
  f  [label="Fetch both strains' genomes\nGFF annotation + protein FASTA from NCBI"];
  same [label="Are A and B the same assembly?", shape=diamond, fillcolor="#fff4d6"];
  idm [label="Identity mapping:\nonly the spelling needs resolving"];
  al [label="Align every protein of A against B, and B against A\n(DIAMOND), keep reciprocal best hits"];
  m  [label="Each other's best hit?", shape=diamond, fillcolor="#fff4d6"];
  y  [label="Same gene in both strains\nPA0001 = PA14_RS00005, 99% identical", fillcolor="#dff0e3"];
  n  [label="No counterpart found:\ngene only in one strain\n(or too diverged / pseudogene)", fillcolor="#f9d9d9"];
  o  [label="Output table: tag in A, tag in B, % identity,\nor 'only in A' / 'only in B'", fillcolor="#dff0e3"];
  s -> w -> ok; ok -> v [label="yes"]; ok -> g [label="no"];
  v -> r; g -> r; r -> f -> same;
  same -> idm [label="yes"]; same -> al [label="no"];
  al -> m; m -> y [label="yes"]; m -> n [label="no"];
  idm -> o; y -> o; n -> o;
}
''')
st.caption("The full method for the *between strains* problem. The same alignment engine is used for species; "
           "see *How species matching works*.")

# ---- live numbers
st.header("What it finds in the real data")
st.write("Computed from the project's tables for every pair of strains of the same species:")
pairs = pd.DataFrame(strains.all_pairs(M))
view = pairs.rename(columns={"species": "Species", "strain_a": "Strain A", "strain_b": "Strain B", "genes_a": "Genes in A",
                             "genes_b": "Genes in B", "shared": "Shared", "only_a": "Only in A", "only_b": "Only in B",
                             "median_identity": "Median % identical", "share_ge_99": "Shared at 99%+"})
st.dataframe(view.style.format({"Shared at 99%+": "{:.0%}", "Median % identical": "{:.1f}"}), hide_index=True, width="stretch")
st.write("Most genes match with almost identical proteins, which is why the hard part here is not the alignment. It is the "
         "genes in the *Only in A / Only in B* columns, and being sure which strain a study used.")

st.subheader("Example: genes that exist only in PA14")
c = strains.compare(M, "paer_pao1", "paer_pa14")
rows = [{"PA14 locus tag": t, "product": M.genes["paer_pa14"].cds[t]["product"]} for t in c["only_b"]
        if "hypothetical" not in M.genes["paer_pa14"].cds[t]["product"]]
hyp = len(c["only_b"]) - len(rows)
st.dataframe(pd.DataFrame(rows[:8]), hide_index=True, width="stretch")
st.caption(f"{len(c['only_b'])} PA14 genes have no PAO1 counterpart; {hyp} of them are annotated only as \"hypothetical protein\". "
           "Several of the named ones are defence systems and immunity factors, typical of the accessory genome "
           "that differs between strains.")

st.subheader("Example: one gene across the strains")
hub = pd.DataFrame(M.hub_table())
row = hub[hub["symbol"] == "dnaA"].iloc[0]
ec = ["ecol_k12", "ecol_39282", "ecol_rel606"]
ex = [("P. aeruginosa PAO1", "paer_pao1", row["paer_pao1"]), ("P. aeruginosa PA14", "paer_pa14", row["paer_pa14"])] + \
     [(G[k]["label"], k, row[k]) for k in ec]
st.dataframe(pd.DataFrame([{"Strain": n, "locus tag of dnaA": t, "older name": ", ".join(M.genes[k].cds[t]["old"])} for n, k, t in ex]),
             hide_index=True, width="stretch")

# ---- caveats
st.header("What to keep in mind")
st.markdown(
    "- **\"No counterpart found\" does not prove the gene is absent.** It may be too diverged for the cut-offs, split into pieces, "
    "or a pseudogene in the other strain. Genes that are not protein-coding are not compared at all.\n"
    "- **The strain of raw-reads-only studies is inferred** from the metadata and cannot be checked against identifiers.\n"
    "- **Two assemblies of the same strain also have different locus tags.** The tool resolves old and current names, but "
    "check which assembly your data used.\n"
    "- **Duplicated genes** can be matched to the wrong copy; use the identity column to spot them.")

st.header("Try it")
nav.link(st, "strains_compare", "Compare two strains: shared genes, strain-specific genes, and an ID translator")
