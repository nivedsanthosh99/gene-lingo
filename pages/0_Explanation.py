"""Streamlit page: beginner's explanation of the bacterial gene mapping, with examples from the project's real data."""
import csv
from pathlib import Path

import pandas as pd
import streamlit as st

from gene_lingo import bacteria
from gene_lingo.ui import hub_frame

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "bacteria"
SURVEY = ROOT / "data" / "osdr_bacterial_survey.tsv"

st.set_page_config(page_title="Explanation - gene-lingo", layout="centered")

if not (DATA / "genes.tsv.gz").exists():
    st.title("Explanation")
    st.error(f"Mapping tables not found in `{DATA}`. Build them with `gene-lingo bacteria build -o data/bacteria`.")
    st.stop()


@st.cache_resource
def get_mapper():
    return bacteria.Mapper(str(DATA))


M = get_mapper()
G, HUB = M.genomes, M.hub
COLOUR = {"Bacillus subtilis": "green", "Pseudomonas aeruginosa": "blue", "Staphylococcus aureus": "orange",
          "Escherichia coli": "violet", "Rhodospirillum rubrum": "red"}


def badge(key, tag):
    """Locus tag as a badge coloured by species (the same colour everywhere on this page)."""
    return f":{COLOUR.get(G[key]['species'], 'gray')}-badge[{tag.replace('_', chr(92) + '_')}]"


hub_df, others, species = hub_frame(M)
core_n = int((hub_df["n_species"] == len(species) + 1).sum())

st.title("How genes are matched across bacterial species")
st.markdown(
    "GeneLab keeps gene-activity measurements from bacteria grown in space and in space-like conditions. To compare "
    "species you first need to know which gene in one species is the same gene in another. This page explains why "
    "that is hard and how the tool solves it. Every example below comes from the project's real data.")
st.markdown(
    "Jump to: [The problem](#the-problem) · [What the data looks like](#what-the-real-data-looks-like) · "
    "[The solution](#the-solution-compare-proteins-not-names) · [Try it](#try-it) · [Limits](#what-to-keep-in-mind)")

# ---- the hook: one gene, six names (live from hub.tsv)
row = hub_df[hub_df["symbol"] == "dnaA"].iloc[0]
shown = [HUB, "paer_pao1", "paer_pa14", "saur_uams1", "ecol_k12", "rrub_11170"]
st.subheader("The same gene, six different names")
st.write("*dnaA* starts DNA copying and is the first gene in each of these genomes. In the data tables it appears as:")
cols = st.columns(3)
for i, k in enumerate(shown):
    tag = row["locus_tag"] if k == HUB else row[k]
    with cols[i % 3]:
        st.caption(G[k]["label"])
        st.markdown(badge(k, tag))

st.header("The goal")
st.write(
    "Each experiment in GeneLab measures how active every gene is, in one species. The aim is to compare experiments "
    "across species, for example to ask which genes respond to microgravity in all five species studied "
    "(*B. subtilis*, *P. aeruginosa*, *S. aureus*, *E. coli* and *R. rubrum*). That needs a table saying which gene "
    "in species A is the same gene as which one in species B.")

st.header("How a gene is written in the data")
st.write("A gene is a stretch of DNA that codes for a protein. In a GeneLab results table each row is one gene and "
         "the first column is its identifier. The same gene can be written four ways:")
st.markdown(
    "| Kind | What it is | Example |\n|---|---|---|\n"
    f"| **Locus tag** | A systematic ID assigned by the genome annotation | {badge(HUB, 'BSU_00010')} |\n"
    "| **Gene symbol** | A short, human-friendly name | `dnaA` |\n"
    f"| **Probe ID** | An ID from an older microarray chip | {badge('paer_pao1', 'PA0001_dnaA_at')} |\n"
    "| **Protein sequence** | The amino-acid string the gene codes for | `MNAYTVSRLALDAGVSV…` |")
st.write("GeneLab tables use the **locus tag**, or the probe ID for arrays.")

# ---- problem
st.header("The problem")
st.write("Look at the six names for *dnaA* above. Three things make matching them hard:")
st.markdown(
    f"1. **Names do not match across species.** No rule turns {badge('ecol_k12', 'b3702')} into {badge(HUB, 'BSU_00010')}.\n"
    f"2. **One genome can have several spellings.** PA14's gene is {badge('paer_pa14', 'PA14_RS00005')} today and was "
    f"{badge('paer_pa14', 'PA14_00010')} in older annotations. *B. subtilis* is written both {badge(HUB, 'BSU_00010')} and "
    f"{badge(HUB, 'BSU00010')}. A database that expects one spelling, such as Ortholuge, misses genes written in another.\n"
    "3. **Symbols cannot rescue you.** The same gene gets different symbols in different species.")

ec = "ecol_k12"
both = agree = 0
examples = []
for _, r in hub_df.iterrows():
    tag = r.get(ec, "")
    if r["symbol"] and tag:
        s2 = M.genes[ec].cds[tag]["symbol"]
        if s2:
            both += 1
            agree += r["symbol"].lower() == s2.lower()
            if r["symbol"].lower() != s2.lower() and r["symbol"] in ("metS", "hprT", "ebfC"):
                examples.append((r["locus_tag"], r["symbol"], tag, s2))
st.write(f"Computed from the project's tables: of {both:,} matched *B. subtilis*–*E. coli* gene pairs that have a symbol "
         f"on both sides, only **{agree / both:.0%}** share the same symbol. Examples:")
st.dataframe(pd.DataFrame(examples, columns=["B. subtilis", "symbol", "E. coli K-12", "symbol "]), hide_index=True, width="stretch")
st.write("On top of that, *P. aeruginosa* is strain PAO1 in some studies and PA14 in others, *E. coli* appears as three "
         "different strains, and some studies have no gene table at all (raw reads only) or only old chip probe IDs.")

# ---- real data
st.header("What the real data looks like")
st.subheader("GeneLab count tables")
st.write("Each row is a gene, each column a sample, each number a count of sequencing reads.")
st.caption("**OSD-554**, *P. aeruginosa* PA14 (first 5 samples)")
st.dataframe(pd.DataFrame({"gene": ["PA14_RS00005", "PA14_RS00010", "PA14_RS00015"], "4.1": [186, 68, 86], "4.2": [184, 270, 89],
                           "4.3": [304, 600, 309], "4.4": [225, 254, 144], "4.6": [413, 343, 239]}), hide_index=True, width="stretch")
st.caption("**OSD-138**, *B. subtilis* 168 (first 3 samples). Highlighted rows are tRNA genes.")
t138 = pd.DataFrame({"gene": ["BSU_24060", "BSU_28509", "BSU_12770", "trnSL-Arg2", "trnS-Ser", "trnY-Glu"],
                     "Bsub_168_cells_FLT_Rep1": [8019, 646, 14, 56, 3852, 1486], "Rep2": [6017, 523, 10, 66, 2569, 1248],
                     "Rep3": [9464, 754, 40, 58, 2135, 968]})
st.dataframe(t138.style.apply(
    lambda r: ["background-color: rgba(255,170,0,.22)" if r["gene"].startswith("trn") else "" for _ in r], axis=1),
    hide_index=True, width="stretch")
st.caption("**OSD-305**, *Serratia liquefaciens* (first 3 samples; imported from GEO, so the first column is called `Locus Tag`)")
st.dataframe(pd.DataFrame({"Locus Tag": ["M495_RS00005", "M495_RS00010", "M495_RS00015"], "Sample 1": [14694, 5451, 9488],
                           "Sample 2": [8601, 3370, 5855], "Sample 3": [5390, 1859, 2942]}), hide_index=True, width="stretch")
st.write("The tRNA rows in OSD-138 show why a protein-based method drops non-coding genes: they have names, not locus "
         "tags, and no protein. The OSD-305 tags follow RefSeq with a prefix of their own (`M495_RS`), so you cannot "
         "guess the species from the format alone.")

st.subheader("Older microarray studies: probe IDs instead of locus tags")
st.write("OSD-14 (*P. aeruginosa* PAO1, Affymetrix chip) lists probe sets. A gene probe contains the PAO1 locus tag. "
         "Control and intergenic probes do not.")
st.dataframe(pd.DataFrame({
    "ProbesetID": ["PA0001_dnaA_at", "ig_901934_901046_at", "AFFX-Bsubtilis_dapB_at"],
    "SYMBOL": ["dnaA", "NA", "NA"], "ENTREZID": ["NA", "NA", "NA"],
    "what it is": ["a gene probe (contains the PAO1 locus tag)", "intergenic region, no gene", "control spike-in probe"]}),
    hide_index=True, width="stretch")
st.warning("The other annotation columns of the OSD-14 and OSD-15 tables are filled with *yeast* gene names (such as "
           "`YEL002C`), the wrong organism. The probe ID is the only usable key.")
st.write("OSD-43 (*E. coli*, Affymetrix) is a mixed chip: one probe set can carry IDs from two *E. coli* strains at once. "
         "OSD-31 and 39 (*R. rubrum*) give a reporter plus its 60-letter oligo sequence.")
st.dataframe(pd.DataFrame({
    "study": ["OSD-43", "OSD-31 / 39"], "ID": ["1759170_s_at", "RruA0001"], "SYMBOL": ["ECs5237|pepA", ""],
    "ENTREZ / extra": ["913804|948791", "protein 83591341, oligo CCATGTCGTTCAAGGAGCAGTTCCGCAGCG…"],
    "what it is": ["one probe set carries IDs of two E. coli strains", "reporter name plus its 60-letter oligo sequence"]}),
    hide_index=True, width="stretch")

st.subheader("The annotation: where locus tags meet proteins")
st.write("For each genome the method downloads an annotation (GFF) and a protein file (FASTA) from NCBI. One PA14 gene, "
         "shortened, shows both spellings of the name side by side, which is how old names are resolved:")
st.code("gene  483  2027  +  ID=gene-PA14_RS00005; gene=dnaA;\n                     locus_tag=PA14_RS00005; old_locus_tag=PA14_00010\n"
        "CDS   …            Parent=gene-PA14_RS00005; protein_id=WP_003109151.1; gene=dnaA", language=None)

# ---- solution
st.header("The solution: compare proteins, not names")
st.write("Related genes in different species code for similar proteins, because they descend from a common ancestor. "
         "A protein's sequence is a fingerprint that does not depend on what anyone called the gene.")

st.subheader("Step 1. Find the right genome for each dataset")
st.write("OSD-554's description names assembly `GCF_000014625.1` (PA14). The tool checks this against the table: 98% of "
         "the IDs in it are locus tags of that assembly. The rest are tRNA and rRNA genes.")

st.subheader("Step 2. Build a protein list labelled with the dataset's own locus tags")
st.write("Old and new spellings both resolve to the same gene, and probe IDs are translated by rule.")
st.markdown(
    "| Written in the data | Resolved to | Rule |\n|---|---|---|\n"
    f"| {badge('paer_pa14', 'PA14_00010')} | {badge('paer_pa14', 'PA14_RS00005')} | old name |\n"
    f"| {badge(HUB, 'BSU00010')} | {badge(HUB, 'BSU_00010')} | underscore removed |\n"
    f"| {badge('paer_pao1', 'PA0001_dnaA_at')} | {badge('paer_pao1', 'PA0001')} | probe ID contains the locus tag |\n"
    f"| {badge('rrub_11170', 'RruA0001')} | {badge('rrub_11170', 'Rru_A0001')} → {badge('rrub_11170', 'RRU_RS00005')} | reporter name, then old name |")

st.subheader("Step 3. Match proteins between species by reciprocal best hit")
st.write("This step answers: *which protein in species B is the counterpart of this protein in species A?* "
         "It has four parts, shown here with one real protein, *dnaA* of *P. aeruginosa* PA14, compared with every "
         "protein of *B. subtilis*.")

EX = ROOT / "data" / "example_alignment.json"
if EX.exists():
    import json
    ex = json.loads(EX.read_text())
    al, q = ex["alignment"], ex["query"]
    fwd, rev = ex["forward"], ex["reverse"]

    st.markdown("**a) A protein is a string of letters.** Each letter is one of 20 amino acids. "
                f"PA14's {q['symbol']} protein is {q['length']} letters long. Two proteins that descend from a common "
                "ancestor keep many of the same letters in the same places, even after hundreds of millions of years.")

    st.markdown("**b) Alignment: line the two proteins up and score the match.** The tool slides one protein against the other, "
                "allowing gaps, to find the arrangement with the best score. This is the real alignment of the two *dnaA* proteins "
                f"(first {al['columns_shown']} of {al['length']} columns; the first line is PA14, the last line *B. subtilis*):")
    legend = "`|` identical letter, `:` similar amino acid, blank = different, `-` = gap"
    lines = []
    for b in al["blocks"]:
        lines += [f"PA14 {b['qpos']:>4}  {b['q']}", f"          {b['m']}", f"BSU  {b['spos']:>4}  {b['s']}", ""]
    st.code("\n".join(lines).rstrip(), language=None)
    st.caption(legend)
    ident = al["identical"] / al["length"]
    st.write(f"Over the whole alignment of {al['length']} columns, {al['identical']} are identical (**{ident:.1%}**, the "
             f"`pident` value), {al['similar']} more are similar, and {al['gaps']} gap characters were needed. "
             "Each pair of letters scores by how often that swap is seen between related proteins (the standard BLOSUM62 table), "
             "so identical or similar letters add to the score and unlikely swaps and gaps subtract. The total is the "
             f"**bit score** ({fwd[0]['bitscore']:.0f}). The **e-value** ({fwd[0]['evalue']:.1e}) is the number of alignments this good "
             "you would expect by pure chance in a search this size, so smaller is better; the tool only keeps hits at 1e-5 or below. "
             "**Coverage** is how much of each protein lines up: a match must cover at least half of both proteins, so a short shared "
             "fragment does not count.")

    st.markdown(f"**c) Best hit: search all {len(M.genes[HUB].cds):,} *B. subtilis* proteins and rank them.** Every protein of PA14 is searched against "
                "every protein of *B. subtilis* with DIAMOND, a program that is much faster than the classic BLAST: it looks for short "
                "exact word matches first and extends only the promising ones (about 3 seconds for a pair of bacterial genomes on 4 "
                f"cores). For {q['tag']} the hits that pass the e-value cut-off, best first:")
    st.dataframe(pd.DataFrame([{"B. subtilis protein": h["subject"], "gene": h["symbol"], "% identical": h["pident"],
                                "bit score": h["bitscore"], "e-value": f"{h['evalue']:.1e}",
                                "covers (PA14 / BSU)": f"{h['query_coverage']:.0%} / {h['subject_coverage']:.0%}",
                                "kept?": "yes, best hit" if i == 0 and h["passes_coverage"] else ("no: covers too little" if not h["passes_coverage"] else "no: not the best")}
                               for i, h in enumerate(fwd)]), hide_index=True, width="stretch")

    st.caption("Coverage here is the share of each protein's letters that line up. The tool's rule divides the alignment length "
               "(which also counts gap columns) by each protein's length; the two differ by a few percent when there are gaps.")
    st.markdown("**d) Reciprocal: ask the question back.** The best hit in one direction is not enough: every protein has a closest relative "
                "in the other species, even when its true counterpart is missing. So the search is repeated from "
                f"*B. subtilis* {ex['reverse_query']} against all PA14 proteins:")
    st.dataframe(pd.DataFrame([{"PA14 protein": h["subject"], "gene": h["symbol"], "% identical": h["pident"], "bit score": h["bitscore"],
                                "e-value": f"{h['evalue']:.1e}", "covers (BSU / PA14)": f"{h['query_coverage']:.0%} / {h['subject_coverage']:.0%}",
                                "kept?": "yes, best hit" if i == 0 and h["passes_coverage"] else "no"} for i, h in enumerate(rev)]),
                 hide_index=True, width="stretch")
    if ex["reciprocal"]:
        st.success(f"{q['tag']} picks {fwd[0]['subject']}, and {fwd[0]['subject']} picks {q['tag']} back. "
                   "Each is the other's best hit, so they are called orthologs. The tiny identity difference between the two "
                   f"directions ({fwd[0]['pident']}% and {rev[0]['pident']}%) is normal: the two searches align slightly differently.")
    with st.expander("Why not simply take the best hit one way?"):
        st.write("Imagine a gene that exists in *P. aeruginosa* but was lost in *B. subtilis*. Searching one way it still gets a "
                 "best hit, a distant relative from a different gene family. Searching back from that relative, its best hit "
                 "is its own true partner, not the first gene, so the pair is rejected. Requiring both directions removes most "
                 "such false matches. It cannot remove all of them: in gene families it can still pick the wrong member "
                 "(see the limits at the bottom).")
else:
    st.write("Compare every protein in species A against species B, and the reverse. Two genes count as the same only if each "
             "is the other's best match.")

st.write("Repeating this for every protein of PA14 gives the ortholog table. Its first rows, live from the project's table "
         "(`pident` is the percentage of identical amino acids):")
pair = M.pair("paer_pa14", HUB)
st.dataframe(pd.DataFrame([(k, v[0], v[1]) for k, v in list(pair.items())[:5]], columns=["PA14", "B. subtilis", "pident"]),
             hide_index=True, width="stretch")
st.write("Even the real *dnaA* pairs are only about 40% identical between these species. That is why exact-name matching "
         "cannot work, and why the tool reports the percentage.")

st.subheader("Step 4. Use *B. subtilis* as a common reference (the hub)")
st.write("Every species is mapped to *B. subtilis*, so all results can be written in its locus tags. Live rows of the hub table:")
view = hub_df[hub_df["symbol"].isin(["dnaA", "rpoB", "recA"])][["locus_tag", "symbol", *others]].rename(
    columns={**{k: G[k]["label"] for k in others}, "locus_tag": G[HUB]["label"]})
st.dataframe(view, hide_index=True, width="stretch")
st.write(f"**{core_n}** *B. subtilis* genes have a counterpart in all five species. This is the **core** gene set, the "
         "genes you can compare across every experiment. Species pairs are also compared directly, which finds genes "
         "the *B. subtilis* route misses.")

st.subheader("Step 5. Make it work for new genomes and new studies")
st.write("**Any genomes:** give an accession or a name such as \"Aliivibrio fischeri ES114\" and the same steps run. "
         "**Any OSDR study:** give an OSD number. The organism and strain come from the study's metadata, the genome is "
         "chosen, and it is checked against the study's table when one exists. Live from the survey of all bacterial "
         "transcription studies in OSDR:")
if SURVEY.exists():
    with open(SURVEY) as fh:
        srv = [r for r in csv.DictReader(fh, delimiter="\t") if r["osd"] in ("OSD-554", "OSD-305", "OSD-356")]
    st.dataframe(pd.DataFrame([{"Study": r["osd"], "Genome chosen": r["genome_label"], "How sure": r["confidence"],
                                "IDs found": r["id_match"] or "n/a (no gene table)"} for r in srv]),
                 hide_index=True, width="stretch")
st.caption("validated = the identifiers in the study's table match the genome; strain = matched by strain name only "
           "(raw reads, nothing to check against).")

# ---- try it
st.header("Try it")
st.write("Type any identifier from the five species' genomes (a current or old locus tag, a probe ID), or pick an example. "
         "The tool finds the genome, resolves the name and lists the same gene in every other genome.")
EXAMPLES = ["PA14_00010", "BSU00010", "PA0001_dnaA_at", "RruA0001", "b3987", "QV15_00005", "ECB_03585", "PA14_RS99999"]
st.session_state.setdefault("ex_input", "PA14_00010")


def pick(x):
    st.session_state["ex_input"] = x


ecols = st.columns(4)
for i, x in enumerate(EXAMPLES):
    ecols[i % 4].button(x, key=f"ex_{x}", on_click=pick, args=(x,), width="stretch")
q = st.text_input("Identifier", key="ex_input").strip()

if q:
    src = M.detect([q])
    if src is None:
        st.warning(f"**{q}** is not a gene in any of the {len(G)} genomes of this project (the tool says so instead of guessing).")
    else:
        tag = M.resolve(q, src)[0]
        how = {"current": "is the current locus tag of", "old": "is an old name for",
               "nounderscore": "is the locus tag written without its underscore for",
               "other": "is an array probe, reporter or Entrez id that resolves to"}[M.how(q, src, tag)]
        st.success(f"{badge(src, q)} {how} {badge(src, tag)}, a gene of **{G[src]['label']}**.")
        res = M.map([q], src, list(G))[0]
        st.markdown(f"***{res['symbol'] or '(no symbol)'}***: {res['product']}. Hub gene (the same gene in *B. subtilis*): "
                    f"{badge(HUB, res['hub_gene']) if res['hub_gene'] else 'none found'}")
        rows = []
        for k in G:
            t = res[k]
            old = ", ".join(M.genes[k].cds[t]["old"]) if t else ""
            pid = None if (k == src or not t) else res[f"{k}_pident"]
            rows.append({"Genome": G[k]["label"] + ("  (your input)" if k == src else ""), "Locus tag": t or "no ortholog found",
                         "Also written as": old, "% identical to your gene": pid})
        st.dataframe(pd.DataFrame(rows).style.format({"% identical to your gene": "{:.1f}"}, na_rep=""),
                     hide_index=True, width="stretch")

# ---- limits
st.header("What to keep in mind")
st.markdown(
    "- **One best match per gene.** In gene families the tool can pick the wrong member. *B. subtilis* `lysS` was matched "
    "to *E. coli* `lysU`. The percent-identity column helps spot these.\n"
    "- **Very diverged genes can be missed.** *R. rubrum*'s `ftsZ` was found in direct comparison but not through the hub.\n"
    "- **Not yet compared with other tools.** The results have not been checked against Ortholuge or OrthoFinder, and the "
    "cut-offs (50% coverage, e-value 1e-5) are conventional rather than tuned.\n"
    "- **Protein-coding genes only.** tRNA, rRNA and other non-coding entries are not mapped.")

st.header("Where to go next")
for col, (path, text, icon) in zip(st.columns(3), [
        ("pages/2_Bacterial_genomes.py", "Bacterial genomes: map your own list", ":material/table_view:"),
        ("pages/3_Any_genomes.py", "Any genomes: your organisms", ":material/biotech:"),
        ("pages/4_Any_OSDR_study.py", "Any OSDR study: by OSD number", ":material/science:")]):
    try:
        col.page_link(path, label=text, icon=icon)
    except st.errors.StreamlitPageNotFoundError:   # page run on its own (tests), not through app.py
        col.write(text)
