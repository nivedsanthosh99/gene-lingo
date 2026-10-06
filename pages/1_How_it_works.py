"""Streamlit page: flowchart + API/format reference."""
import pandas as pd
import streamlit as st

from gene_lingo import ID_TYPES
from gene_lingo.species import SPECIES

st.set_page_config(page_title="How it works - gene-lingo", layout="wide")
st.title("How it works")

FLOW = r'''
digraph G {
  rankdir=TB; bgcolor="transparent"; nodesep=0.35; ranksep=0.45;
  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=11, fillcolor="#eef2f7", color="#6b7c93"];
  edge [fontname="Helvetica", fontsize=10, color="#6b7c93"];

  in   [label="User input\n(text box or CSV/TSV upload, first column)\n+ input species, output species,\ninput ID type, output ID type", fillcolor="#dff0e3"];
  norm [label="Clean list\n(strip, drop blanks, de-duplicate)"];
  res  [label="Species resolver (species.py)\ncommon name / Latin / Ensembl name / taxid\n-> (NCBI taxid, Ensembl name)"];
  same [label="Same species?", shape=diamond, fillcolor="#fff4d6"];

  conv1 [label="convert()\nid -> output ID type", fillcolor="#e3e8fb"];

  oin  [label="Input type usable for orthology?\n(auto, symbol, ensembl, entrez, hgnc)", shape=diamond, fillcolor="#fff4d6"];
  err  [label="Error:\nUniProt/RefSeq/name input\nneeds same-species conversion", fillcolor="#f9d9d9"];
  toens[label="convert()\ninput -> Ensembl gene ID\n(source species)", fillcolor="#e3e8fb"];
  ens  [label="Ensembl Compara homology\nper Ensembl gene ID", fillcolor="#e3e8fb"];
  filt [label="Filter by orthology type\n(all / one2one / one2many / many2many)"];
  conv2[label="convert()\northolog Ensembl ID -> output ID type\n(target species; skipped if output = ensembl)", fillcolor="#e3e8fb"];

  mg   [label="MyGene.info\nPOST https://mygene.info/v3/query\nform: q, scopes, fields, species=taxid, size=5\nbatches of 500", shape=cylinder, fillcolor="#fde2c8"];
  ensapi[label="Ensembl REST\nGET https://rest.ensembl.org/homology/id/{species}/{ENSG}\n?type=orthologues&target_species=...&sequence=none", shape=cylinder, fillcolor="#fde2c8"];

  out  [label="Rows: input, input_species, output_species,\noutput_id, orthology_type, perc_id\n(blank output_id = not mapped)", fillcolor="#dff0e3"];
  ui   [label="Results table (unmapped rows highlighted)\n+ 'N of M mapped' + Download CSV", fillcolor="#dff0e3"];

  in -> norm -> res -> same;
  same -> conv1 [label="yes"];
  same -> oin   [label="no"];
  oin -> err    [label="no"];
  oin -> toens  [label="yes"];
  toens -> ens -> filt -> conv2 -> out;
  conv1 -> out -> ui;

  conv1 -> mg [style=dashed, dir=both];
  toens -> mg [style=dashed, dir=both];
  conv2 -> mg [style=dashed, dir=both];
  ens -> ensapi [style=dashed, dir=both];
}
'''
st.graphviz_chart(FLOW, width="stretch")
st.caption("Solid arrows: data flow. Dashed arrows: calls to external APIs. Yellow diamonds: decisions.")

st.header("External APIs")
st.dataframe(pd.DataFrame([
    {"API": "MyGene.info", "Used for": "ID conversion (all ID types, alias/old symbols)",
     "Request": "POST /v3/query (form-encoded): q=<comma-joined ids>, scopes, fields, species=<taxid>, size=5",
     "Response": "JSON list of hits; each has `query`, requested fields, or `notfound: true`",
     "Notes": "500 ids per request; no API key"},
    {"API": "Ensembl REST", "Used for": "Orthologs (Compara)",
     "Request": "GET /homology/id/{source ensembl species}/{Ensembl gene id}?type=orthologues&target_species=...&sequence=none",
     "Response": "JSON `data[].homologies[]`: `type`, `target.id`, `target.species`, `target.perc_id`, `target.perc_pos`",
     "Notes": "HTTP 400 treated as 'no result'; HTTP 429 retried using Retry-After; no API key"},
]), hide_index=True, width="stretch")

st.header("ID types (MyGene field used)")
st.dataframe(pd.DataFrame(
    [{"ID type (UI)": k, "MyGene field": v} for k, v in ID_TYPES.items()]
    + [{"ID type (UI)": "auto (input only)",
        "MyGene field": "tries scopes: symbol, alias, entrezgene, ensembl.gene, HGNC, uniprot, refseq, name"}]),
    hide_index=True, width="stretch")

st.header("Species formats")
st.write("Any of these are accepted as input. Output rows currently use the **Ensembl name**.")
st.dataframe(pd.DataFrame(
    [{"Common name": c, "NCBI taxid": t, "Ensembl name": e, "Latin name": e.replace("_", " ").capitalize()}
     for c, t, e in SPECIES]), hide_index=True, width="stretch")

st.header("Output columns")
st.dataframe(pd.DataFrame([
    {"Column": "input", "Meaning": "Gene exactly as entered"},
    {"Column": "input_species / output_species", "Meaning": "Ensembl species names"},
    {"Column": "output_id", "Meaning": "Gene in the chosen output ID type; blank if not mapped"},
    {"Column": "orthology_type", "Meaning": "Ensembl type, e.g. ortholog_one2one (blank for same-species conversion)"},
    {"Column": "perc_id", "Meaning": "Percent identity of the target to the query, from Ensembl"},
]), hide_index=True, width="stretch")
