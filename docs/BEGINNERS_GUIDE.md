# Mapping bacterial genes across species: a beginner's guide

## The goal

NASA's GeneLab/OSDR holds experiments where bacteria were grown in spaceflight or simulated microgravity. Each
experiment measures how active every gene is, in one species. The aim is to compare experiments across species,
for example to ask which genes respond to microgravity in all five species. That requires knowing which gene in
species A is "the same gene" as one in species B.

## How a gene is written in the data

A gene is a stretch of DNA that codes for a protein. In a GeneLab results table each row is one gene and the first
column is its identifier. A gene can be written in several ways:

| Kind | What it is | Example |
|---|---|---|
| **Locus tag** | A systematic ID assigned by the genome annotation | `BSU_00010` |
| **Gene symbol** | A short human-friendly name | `dnaA` |
| **Probe ID** | An ID from an older microarray chip | `PA0001_dnaA_at` |
| **Protein sequence** | The amino-acid string the gene codes for | `MNAYTVSRLALDAGVSV...` |

GeneLab tables use the **locus tag** (or, for arrays, the probe ID).

## The problem

The same gene has a different name in every species. The gene `dnaA` (it starts DNA copying) is the first gene in
all of these genomes:

| Species | Locus tag in the data |
|---|---|
| *B. subtilis* | `BSU_00010` |
| *P. aeruginosa* PAO1 | `PA0001` |
| *P. aeruginosa* PA14 | `PA14_RS00005` |
| *S. aureus* | `QV15_00005` |
| *E. coli* K-12 | `b3702` |
| *R. rubrum* | `RRU_RS00005` |

Three things make this hard:

1. **Names don't match across species.** No rule turns `b3702` into `BSU_00010`.
2. **One genome can have several spellings.** PA14's gene is `PA14_RS00005` today but was `PA14_00010` in older
   annotations; *B. subtilis* is written both `BSU_00010` and `BSU00010`. A database expecting one scheme (such as
   Ortholuge) misses genes written in another.
3. **Symbols can't rescue you.** The same gene is `metS` in *B. subtilis* but `metG` in *E. coli*, and `hprT` vs `hpt`.
   Only 49 % of matched *B. subtilis*–*E. coli* pairs even share a symbol.

Further complications: *P. aeruginosa* is strain PAO1 in some studies and PA14 in others; *E. coli* appears as
three strains; some studies have no gene table at all (raw reads only) or old chip probe IDs.

## What the real data looks like

### The input: GeneLab tables

Each row is a gene, each column a sample, each number a count of sequencing reads (OSD-554, *P. aeruginosa* PA14,
first 5 columns):

```
"","4.1","4.2","4.3","4.4","4.6"
"PA14_RS00005",186,184,304,225,413
"PA14_RS00010",68,270,600,254,343
"PA14_RS00015",86,89,309,144,239
```

The same layout for a different organism and annotation (OSD-138, *B. subtilis* 168). The last rows show why
non-coding genes drop out of a protein-based method: they have names, not locus tags, and no protein.

```
"","Bsub_168_cells_FLT_Rep1","Bsub_168_cells_FLT_Rep2","Bsub_168_cells_FLT_Rep3"
"BSU_24060",8019,6017,9464
"BSU_28509",646,523,754
...
"trnSL-Arg2",56,66,...
"trnS-Ser",3852,2569,...
```

A study imported from GEO (OSD-305, *Serratia liquefaciens*) uses a header called `Locus Tag`; the tags follow
RefSeq, with a prefix of its own (`M495_RS`), so you cannot guess the species from the format alone:

```
Locus Tag      Sample 1  Sample 2  Sample 3
M495_RS00005   14694     8601      5390
M495_RS00010   5451      3370      1859
```

### Older microarray studies: probe IDs, not locus tags

OSD-14 (*P. aeruginosa* PAO1, Affymetrix chip) lists probe sets. `PA0001_dnaA_at` contains the PAO1 locus tag; control
probes and intergenic probes do not:

```
ProbesetID                SYMBOL   ENTREZID
PA0001_dnaA_at            dnaA     NA
ig_901934_901046_at       NA       NA        <- intergenic region, no gene
AFFX-Bsubtilis_dapB_at    NA       NA        <- control spike-in probe
```

Beware: the other annotation columns of the OSD-14/15 tables are filled with *yeast* gene names (`YEL002C`), the
wrong organism, so the probe ID is the only usable key.

OSD-43 (*E. coli*, Affymetrix) is a mixed chip: one probe set can carry IDs from two *E. coli* strains at once:

```
ProbesetID     SYMBOL         ENTREZID
1759170_s_at   ECs5237|pepA   913804|948791
```

OSD-31/39 (*R. rubrum*) give reporters plus a 60-letter oligo sequence:

```
RruA0001   83591341   CCATGTCGTTCAAGGAGCAGTTCCGCAGCG...
```

### The annotation: where locus tags meet proteins

The method downloads, for each genome, an annotation (GFF) and a protein file (FASTA) from NCBI. One gene in the
PA14 annotation (shortened) shows both spellings of the name side by side, which is how old names are resolved:

```
gene  483  2027  +  ID=gene-PA14_RS00005; gene=dnaA; locus_tag=PA14_RS00005; old_locus_tag=PA14_00010
CDS   ...  Parent=gene-PA14_RS00005; protein_id=WP_003109151.1; gene=dnaA; ...
```

The protein FASTA then holds the amino-acid sequence for `WP_003109151.1`, which is attached to `PA14_RS00005`.

## The solution: compare proteins, not names

Related genes in different species code for similar proteins because they descend from a common ancestor. A
protein's amino-acid sequence is a fingerprint that does not depend on what anyone called the gene.

**1. Find the right genome for each dataset.** Example: the OSD-554 description names assembly `GCF_000014625.1`
(PA14), and 98 % of the IDs in its table are locus tags of that assembly (the rest are tRNA/rRNA genes).

**2. Build a list of proteins, each labelled with the dataset's locus tag.** Old and new spellings both resolve
(`PA14_00010` and `PA14_RS00005` → the same gene). Probe IDs are translated by rule: `PA0001_dnaA_at` → `PA0001`,
`RruA0001` → `Rru_A0001`.

**3. Match proteins between species by reciprocal best hit.** Compare every protein in species A against species
B and the reverse; two genes are "the same" only if each is the other's best match. Output for PA14 vs *B. subtilis*
(`pident` = % identical amino acids):

```
paer_pa14      bsub_168   pident
PA14_RS00005   BSU_00010  38.4
PA14_RS00010   BSU_00020  26.5
PA14_RS00020   BSU_00060  48.1
```

Even the real `dnaA` pairs are only about 40 % identical between these species, which is why exact-name matching
cannot work and why the method uses sequence similarity and shows the percentage.

**4. Use *B. subtilis* as a common reference (hub).** Every species is mapped to *B. subtilis*, so all results can be
written in its locus tags. One row of the result table (`hub.tsv`) lists a gene across all genomes:

```
BSU_00010  dnaA  chromosomal replication initiator informational ATPase
  PAO1: PA0001   PA14: PA14_RS00005   S. aureus: QV15_00005   E. coli K-12: b3702
  E. coli NCTC 86: E4V11_RS00005   E. coli REL606: ECB_RS18960   R. rubrum: RRU_RS00005
BSU_01070  rpoB  RNA polymerase (beta subunit)
  PAO1: PA4270   PA14: PA14_RS03545   S. aureus: QV15_02465   E. coli K-12: b3987   R. rubrum: RRU_RS13915
```

488 *B. subtilis* genes have a counterpart in all five species. This is the **core** gene set, the genes you can
compare across every experiment. Species pairs are also compared directly, which finds genes the *B. subtilis*
route misses.

**5. Make it work for new data.**
- *Any genomes:* give an accession or a name such as "Aliivibrio fischeri ES114" and the same steps run.
- *Any OSDR study:* give an OSD number; the organism and strain come from the metadata and the genome is chosen,
  and checked against the study's table when one exists. For two studies this gave (from the survey file):

```
OSD-554  P. aeruginosa UCBPP-PA14   validated: 2931 of 3000 sampled IDs found in the genome
OSD-305  S. liquefaciens ATCC 27592 validated: 2890 of 3000 sampled IDs found in the genome
OSD-356  A. fischeri ES114          strain name matched (raw reads only, nothing to validate against)
```

### Using the result

Translating the three genes below from OSD-554 (PA14) into *B. subtilis* terms: an old-style name, a current name,
and one that does not exist.

```
$ gene-lingo bacteria map --from OSD-554 --to OSD-138 -d data/bacteria
input          locus_tag      B. subtilis 168   pident
PA14_00010     PA14_RS00005   BSU_00010         38.4      <- old name resolved, then mapped
PA14_RS23350   PA14_RS23350   BSU_15290         49.5      (ftsZ)
PA14_RS99999                                              <- not a gene in this genome
```

The same gene to *E. coli*: `PA14_RS00005` → `E4V11_RS00005` (64.1 % identical).

I surveyed all 28 bacterial transcription studies in OSDR this way. For the 19 with no gene table, the genome
choice cannot be verified.

## What to keep in mind

- The method gives one best match per gene, so in gene families it can pick the wrong member. Example: *B.
  subtilis* `lysS` was matched to *E. coli* `lysU`. The percent-identity column helps spot these.
- Very diverged genes can be missed: *R. rubrum*'s `ftsZ` was found in direct comparison but not through the hub.
- The results have not been compared against other tools such as Ortholuge or OrthoFinder, and the cut-offs
  (50 % coverage, e-value 1e-5) are conventional rather than tuned.
- Only protein-coding genes are mapped.

More detail: `docs/METHOD.md` (technical method) and `README.md` (how to run it).
