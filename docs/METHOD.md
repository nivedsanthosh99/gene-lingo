# Mapping bacterial gene identifiers across OSDR datasets: method

## Summary

Gene identifiers in the OSDR bacterial datasets cannot be matched by name: the same genome is annotated with
several locus-tag schemes, strains differ between studies, and gene symbols are inconsistent between species.
The method therefore never compares names. It (1) establishes which genome each dataset's identifiers belong to,
(2) builds a proteome from that genome's annotation, so every protein carries the locus tag the dataset uses,
(3) finds orthologs between proteomes by reciprocal best protein hits, and (4) expresses all results through a
hub genome (*B. subtilis* 168) as well as directly for every pair of genomes. Results are therefore keyed on the
identifiers that appear in the GeneLab count tables. The same procedure is available for arbitrary genomes and,
with the genome inferred automatically, for any OSDR study.

## 1. Why name-based mapping fails

Three things were established from the OSDR files themselves:

- **One genome, several locus-tag schemes.** Examples: *P. aeruginosa* PA14 `PA14_RS00005` (current RefSeq) vs
  `PA14_00010` (original), and PAO1 `PA0001`; *B. subtilis* `BSU_00010` vs `BSU00010`; *R. rubrum* `RRU_RS00005` vs
  `Rru_A0001`. A database keyed to one scheme (e.g. Ortholuge) silently misses genes written in another.
- **Strains differ between and within species.** *P. aeruginosa* is PAO1 in OSD-14/15 but PA14 in OSD-554; the
  *E. coli* studies are K-12 MG1655 (OSD-43), ATCC 4157 / NCTC 86 (OSD-95) and B REL606 (OSD-728).
- **Gene symbols are not a usable key.** Among reciprocal best hits between *B. subtilis* and *E. coli* K-12 only
  49 % of pairs with a symbol on both sides have the same symbol (e.g. `metS`/`metG`, `hprT`/`hpt`).

Several OSDR bacterial studies also have no processed gene table at all (OSD-215, 216, 728, 31, 39 and most others),
and the array studies report probe ids rather than locus tags (OSD-14/15 Affymetrix `PA0001_dnaA_at`; OSD-43
Affymetrix probe sets with Entrez ids; OSD-31/39 reporters `RruA0001`).

## 2. Step 1: tie each dataset to a genome

**Curated set (five species, eight genomes).** For each study the assembly was taken from the OSDR processing
description (e.g. OSD-554: `GCF_000014625.1_ASM1462v1`; OSD-138/185: Ensembl `ASM904v1`, the same assembly as
RefSeq `GCF_000009045.1`; OSD-145: `GCA_000788115.1`; OSD-95: `GCF_900699165.1`) and then *checked* against the
identifiers in the study's count table: 98 % (OSD-554), 98 % (138, 185), 94 % (145) and 92 % (95) of the table's
identifiers are locus tags of that assembly. The remainder are non-coding entries (tRNA, rRNA, pseudogenes),
which have no protein and are outside this method. For studies without a table (OSD-215/216 *B. subtilis*,
OSD-728 *E. coli* B REL606, OSD-31/39 *R. rubrum*) the genome follows from the strain in the metadata.

**Any OSDR study (automatic; section 6).** The same decision is made from the OSDR metadata and, where one exists,
the processed table.

## 3. Step 2: proteomes keyed by dataset locus tag

For each genome the NCBI genomic GFF3 and protein FASTA are downloaded. A gene is kept if its CDS record has both a
`locus_tag` and a `protein_id` present in the FASTA; the protein sequence is attached to that locus tag (one protein
per locus tag, first CDS record). Pseudogenes and RNA genes therefore drop out. Each gene also keeps its symbol,
product, `old_locus_tag` values and Entrez GeneID (from `Dbxref`). Resulting counts: PAO1 5,572; PA14 5,891;
*B. subtilis* 168 4,237; *S. aureus* UAMS-1 2,653; *E. coli* K-12 4,290, NCTC 86 4,529, REL606 4,204;
*R. rubrum* 3,870.

An alias table maps every identifier form to the canonical (current) locus tag: the current tag, each
`old_locus_tag`, and the tag with underscores removed. `PA14_00010`, `PA14_RS00005`, `BSU00010` and `BSU_00010`
therefore all resolve. Probe and reporter ids are translated by rule, each verified against data:

| Study | Id form | Rule | Check |
|---|---|---|---|
| OSD-14, 15 | `PA0001_dnaA_at`, `PA0005_at` | leading `PA####` is the PAO1 locus tag; `ig_…` and `Pae_…` control probes have none | 98 % of the table's 5,900 probe sets resolve |
| OSD-31, 39 | `RruA0001` | → `Rru_A0001` (an `old_locus_tag`) | the 60-mer oligo of 3,786 reporters was located in the genome: 96.6 % fall inside the CDS of the assigned locus tag |
| OSD-43 | Entrez id(s), `id1\|id2` | GeneID → locus tag via the K-12 GFF | K-12 genes only; O157/CFT073-specific probes do not map |

## 4. Step 3: orthologs by reciprocal best hits

For every pair of proteomes A and B:

1. DIAMOND 2.1.9 `blastp`, `--very-sensitive`, e-value ≤ 1e-5, up to 5 hits per query, run in **both** directions.
2. Keep only alignments whose length covers **≥ 50 % of the query and ≥ 50 % of the subject** (alignment length /
   protein length; the alignment length also counts gap columns, so it can slightly exceed the number of aligned residues). No identity cut-off is applied; percent identity is carried through as a quality column.
3. For each query take the highest-bitscore surviving hit.
4. A pair (a, b) is an ortholog if b is a's best hit **and** a is b's best hit.

This yields one-to-one pairs by construction. All 28 pairs of the eight curated genomes are computed. Results are
cached by assembly accession and thresholds, so adding a genome only computes the new pairs.

**Hub and core.** *B. subtilis* 168 is the hub: `hub.tsv` has one row per *B. subtilis* gene and the ortholog in each
other genome. Each other genome has an ortholog for 952 (*R. rubrum*) to 1,391 (*S. aureus*) of the 4,237 hub genes.
The **core** is the set of hub genes with an ortholog in at least one genome of every species, i.e. tolerant of a
missing strain: 488 genes for the five species (481 if all eight genomes are required). Any pair can also be
derived *through* the hub (a → hub → b, % identity = the lower of the two); in `--hub-only` builds, which need only
n−1 comparisons, non-hub pairs are always derived this way.

## 5. Why the hub is not the only route

Direct and hub-derived mappings were compared for all 21 non-hub pairs, over genes that have a hub ortholog on
both sides. They give the same gene for 94.7 % overall (≈ 89–90 % for pairs with *S. aureus*, 93–95 % for the
other cross-species pairs, 99.6–100 % between strains of one species). Direct pairs
are the better mapping and cover genes absent from *B. subtilis*: *P. aeruginosa* PA14 has 1,153 hub orthologs but 1,918
direct orthologs with *E. coli* K-12. The hub is for a consistent coordinate system and the core, not a substitute
for direct comparison.

## 6. Any genome, and any OSDR study

**Any genomes.** A genome is an NCBI assembly accession, an organism/strain name resolved against a slim local copy
of NCBI's RefSeq assembly summary (current bacterial/archaeal assemblies that are complete or
reference/representative; ranked by exact strain name, then reference > representative > complete > newest), a
`taxid:`, or user-supplied GFF + protein FASTA. Steps 2–4 then apply unchanged.

**Any OSDR study.** For each study the metadata gives organism, NCBI taxon and strain; any assembly accession
mentioned in the processing text (also from an Ensembl annotation file name) is extracted. Candidate genomes are
ordered: (1) an assembly named in the metadata; (2) organism + strain by name, including a parent strain named in
the strain text (e.g. "chi3339 (mouse-passaged derivative of SL1344)" → SL1344); (3) the study's taxon and the species'
reference genome; (4) an arbitrary strain of the taxon as a last resort. If the study has a processed table
(unnormalised counts, a GEO raw-count table, or a microarray table), up to 3,000 identifiers are scored against the
first six candidates (current/old locus tags; for arrays the probe id with the `_at` suffix and trailing annotation
stripped, and Entrez ids). The first candidate whose match rate is at least 45 % **and** within 15 points of the best
is chosen (**validated**). The margin matters for multi-strain arrays: the OSD-43 chip recognises *E. coli* O157:H7
Sakai (78 %) better than K-12 (68 %), but the study is K-12 and the better-named candidate wins. Without a table
the best name-based candidate is returned, graded *stated* (named in the metadata), *strain*, *species* or *taxon*.

Applied to the 28 bacterial transcription studies in OSDR (17 distinct genomes): 9 validated, 1 stated, 14
matched by strain name, 4 species-level. The 19 that are not validated are studies without a processed gene table, so the
choice cannot be verified.

## 7. Evidence that it works

- **Same-species consistency.** PAO1 ↔ PA14 gives 5,346 reciprocal hits (median 99.6 % identity); the three
  *E. coli* strains give 3,766–3,889 (median 100 %).
- **Hub vs direct agreement:** 94.7 % (section 5).
- **Probe rules:** 96.6 % of *R. rubrum* oligos lie in the assigned gene (section 3).
- **Known genes:** rpoB, rpoC, gyrA, gyrB, recA, dnaK, groEL, rpsC, rplB, tufA, atpA, dnaN, sigA and sdhA are in the
  core; the sporulation regulator spo0A is not.
- **Generality:** the automatic OSDR pipeline reproduces the curated *B. subtilis*–PA14 ortholog set exactly
  (1,153 of 1,153 pairs).
- **Tests:** 23 offline unit/UI tests cover parsing, probe rules, ranking, genome choice, caching and hub
  derivation (the last found a real bug that full builds had masked).

## 8. Limitations

- **Not benchmarked against other ortholog resources.** Ortholuge and OrthoFinder were not run, so agreement with
  them is unknown (Ortholuge's site was not reachable from the build environment). The thresholds (50 %
  coverage, e-value 1e-5, no identity filter) are conventional, not tuned.
- **One ortholog per gene.** Reciprocal best hit can pick the wrong member of a paralog family, e.g. *B. subtilis*
  `lysS` → *E. coli* `lysU` (K-12 has both), `clpC` → `clpB`. Use the percent-identity column to triage.
- **The 50 % coverage rule drops divergent or multi-domain genes.** *R. rubrum* `ftsZ` has no hub hit although the
  direct *P. aeruginosa*–*R. rubrum* pair finds it; hub tables and the core are conservative.
- **Hub bias.** Genes absent from *B. subtilis* are not in the hub table or core; the core size (488 for these five
  species) depends on the hub and the thresholds.
- **Protein-coding genes only.** tRNA, rRNA and other non-coding entries (e.g. OSD-138's `trn…`, `rrn…` names) are not
  mapped.
- **Genomes for raw-only studies are inferred.** A different sub-strain changes strain-specific genes only, but
  the genome is not verified without a table (OSD-11, 31, 323, 810 are species-level). Even for validated
  studies 2–8 % of table identifiers are not recognised; these are mostly non-coding entries, and annotation drift
  between the GTF a study used and today's NCBI GFF may account for some.
- **Array coverage.** Probe rules exist for the PAO1 chip and the *R. rubrum* array and Entrez ids for OSD-43;
  other arrays (e.g. the *Salmonella* chip, OSD-11) would need a rule.

## 9. Reproduce

```
sudo apt install diamond-aligner
gene-lingo bacteria build -o data/bacteria        # curated eight genomes, ~1.5 min
gene-lingo osdr survey -o survey.tsv              # genome choice for every bacterial OSDR study
gene-lingo osdr build OSD-305 OSD-356 OSD-554 -o out/
gene-lingo genomes build "Bacillus subtilis 168" GCF_000011805.1 -o out/
gene-lingo bacteria map -d out/ --from OSD-305 --to OSD-554 M495_RS00005
```

Code: `gene_lingo/bacteria.py` (proteomes, reciprocal best hits, hub, mapping), `genomes.py` (assembly lookup),
`osdr.py` (OSDR metadata and genome choice). Outputs: `data/bacteria/` (curated), `data/osdr_all/` (all studies),
`data/osdr_bacterial_survey.tsv`.
