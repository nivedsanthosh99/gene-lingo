"""Streamlit page: ortholog / locus-tag mapping for any bacterial genomes, not only the OSDR ones."""
import streamlit as st

from gene_lingo import genomes
from gene_lingo.ui import METHOD, build_and_show

st.set_page_config(page_title="Any genomes - gene-lingo", layout="wide")
st.title("Any genomes")
st.caption("Map genes between any bacterial or archaeal genomes. Give organism names, NCBI assembly accessions or your "
           "own annotation files; genes are matched by protein sequence (reciprocal best hits), so any locus-tag "
           "convention works.")

SPECS = st.text_area(
    "Genomes, one per line", height=150, key="any_specs",
    value="Bacillus subtilis 168\nEscherichia coli K-12 MG1655\nPseudomonas aeruginosa PAO1",
    help="Accepted: an assembly accession (GCF_000009045.1, GCA_...), a name with optional strain "
         "(Aliivibrio fischeri ES114), taxid:1280, or local files annotation.gff.gz,proteins.faa.gz,Label")


@st.cache_resource(show_spinner="Loading the NCBI assembly index (downloaded once, ~250 MB)...")
def index():
    return genomes.assembly_index()


registry, errors = {}, []
if SPECS.strip():
    idx = index()
    for line in SPECS.splitlines():
        if line.strip():
            try:
                k, spec = genomes.parse_genome_arg(line, index=idx)
                registry.setdefault(k, spec)
            except (ValueError, OSError) as e:
                errors.append(f"`{line.strip()}`: {e}")
for e in errors:
    st.error(e)
if registry:
    st.dataframe([{"Genome": g["label"], "Species": g["species"], "Assembly": g["acc"] or "local files"}
                  for g in registry.values()], hide_index=True, width="stretch")
    c1, c2 = st.columns(2)
    hub = c1.selectbox("Hub genome", list(registry), format_func=lambda k: registry[k]["label"],
                       help="Everything is also expressed in this genome's locus tags; the core is defined on it.")
    hub_only = c2.checkbox("Only compare each genome with the hub (faster for many genomes)", value=len(registry) > 6)
    build_and_show(registry, hub, hub_only, "any", METHOD)
