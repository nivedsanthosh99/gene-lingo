"""Streamlit page: give any OSDR study ids; find each study's genome and map genes between them."""
import re

import streamlit as st

from gene_lingo import genomes, osdr
from gene_lingo.ui import METHOD, build_and_show

st.set_page_config(page_title="Any OSDR study - gene-lingo", layout="wide")
st.title("Any OSDR study")
st.caption("Enter OSDR study numbers. For each, the organism and strain are read from the OSDR metadata, the genome "
           "is chosen (and, when the study has a processed gene table, checked against the identifiers in it), "
           "then genes are mapped between the studies' genomes.")

txt = st.text_input("OSDR studies", "OSD-554, OSD-138, OSD-305", help="Numbers or OSD-/GLDS- ids, separated by commas or spaces")
override = st.text_area("Override genomes (optional), one per line: `OSD-305 = GCF_000422085.1` or `OSD-90 = Mycobacterium marinum M`", height=70)
studies = [osdr.osd_number(x) for x in re.split(r"[\s,;]+", txt) if x]

if st.button("Look up studies", type="primary", disabled=not studies):
    try:
        with st.status("Reading OSDR metadata and checking genomes...", expanded=True) as status:
            idx = genomes.assembly_index()
            registry, picks = osdr.genomes_for(studies, progress=status.write, index=idx)
            status.update(label="Done", state="complete", expanded=False)
        over = {}
        for line in override.splitlines():
            if "=" in line:
                o, spec = line.split("=", 1)
                over[osdr.osd_number(o)] = genomes.parse_genome_arg(spec, index=idx)
        for n, (key, spec) in over.items():
            if n in picks:
                registry[picks[n]["acc"]]["osds"].remove(n)
                picks[n].update(acc=key, label=spec["label"], species=spec["species"], confidence="you chose it", evidence="override")
                registry.setdefault(key, dict(spec, osds=[]))["osds"].append(n)
        registry = {k: v for k, v in registry.items() if v["osds"]}
        st.session_state["osdr_lookup"] = (studies, registry, picks)
    except Exception as e:  # network / unknown study / no genome
        st.error(str(e))

state = st.session_state.get("osdr_lookup")
if state and state[0] == studies:
    _, registry, picks = state
    st.dataframe([{"Study": f"OSD-{n}", "Organism": "; ".join(p["info"]["organisms"]), "Strain": "; ".join(p["info"]["strains"]),
                   "Gene table": p["info"]["tables"][0][0] if p["info"]["tables"] else "raw reads only",
                   "Genome": f"{p['label']} ({p['acc']})", "Confidence": p["confidence"], "Why": p["evidence"]}
                  for n, p in picks.items()], hide_index=True, width="stretch")
    st.caption("Confidence: validated = identifiers in the study's table match the genome; stated = assembly named in the "
               "study metadata; strain / species / taxon = best match by name, not checked against identifiers (raw-reads-only studies).")
    if len(registry) < 2:
        st.info("These studies share one genome; add a study of another organism or strain to map between them.")
    else:
        hub = next((k for k, g in registry.items() if g["species"] == "Bacillus subtilis"), next(iter(registry)))
        c1, c2 = st.columns(2)
        hub = c1.selectbox("Hub genome", list(registry), index=list(registry).index(hub), format_func=lambda k: registry[k]["label"])
        hub_only = c2.checkbox("Only compare each genome with the hub (faster)", value=len(registry) > 6, key="osdr_hubonly")
        build_and_show(registry, hub, hub_only, "osdr_any", METHOD)
