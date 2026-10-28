# gene-lingo

Two things, both via free public APIs (no keys):

1. **Convert gene names/IDs** between nomenclatures (symbol, Entrez, Ensembl, HGNC, UniProt, RefSeq, full name) using [MyGene.info](https://mygene.info). Old/alias symbols resolve too; input type is auto-detected.
2. **Find orthologs** between species using Ensembl Compara ([REST](https://rest.ensembl.org)), with symbols for the results.

```python
from gene_lingo import convert, orthologs
convert(["TP53", "BRCA1"], to="ensembl")               # {'TP53': ['ENSG00000141510'], ...}
convert(["7157"], to=["symbol", "uniprot"], from_type="entrez")
orthologs(["TP53"], source="human", target="mouse")    # [{'id': 'ENSMUSG...', 'symbol': 'Trp53', 'type': 'ortholog_one2one', ...}]
```

CLI (TSV output; genes as args, `-f file`, or stdin):

```
pip install -e .
gene-lingo convert TP53 BRCA1 --to entrez
gene-lingo orthologs TP53 -s human -t zebrafish --type one2one
```

Species: human, mouse, rat, zebrafish, fly, worm, yeast, chicken, pig, cow, dog, macaque, chimpanzee, frog (or Latin name / taxid); extend `gene_lingo/species.py`.

Tests: `pip install pytest && pytest` (network is mocked).

## Web UI

```
pip install -e ".[ui]"
streamlit run app.py
```

The sidebar has a **How it works** page with the full flowchart, API calls and formats, and three bacterial pages (**Bacterial genomes**, **Any genomes**, **Any OSDR study**, described below). On Streamlit Community Cloud `packages.txt` installs `diamond-aligner` for the build step.

Pick input/output species and ID types, paste genes or upload a CSV, and download the results table as CSV. Same species = ID conversion; different species = ortholog lookup, with the orthologs shown in your chosen output ID type.

## Bacterial locus-tag mapping

OSDR bacterial studies use locus tags from different annotations (`PA14_RS00005` vs old `PA14_00010` vs PAO1 `PA0001`, `BSU_00010` vs `BSU00010`, `RRU_RS00005` vs `Rru_A0001`, ...), so name-based lookups, and tools keyed to one annotation, miss many genes. `gene_lingo.bacteria` instead

1. ties each study to the genome its identifiers come from (NCBI GFF + protein FASTA),
2. computes reciprocal best DIAMOND hits (both proteins ≥50 % covered, e-value ≤1e-5) between every pair of genomes, and against a hub genome (default *B. subtilis* 168) so everything can also be expressed in one genome's locus tags and a core gene set defined,
3. resolves old/alternative locus-tag spellings and array probe ids to the canonical tag.

Three ways in, from most to least curated. All need `diamond` only to *build* tables (`sudo apt install diamond-aligner`); mapping and the web UI read the finished tables.

### 1. The curated OSDR set (`gene-lingo bacteria`, page *Bacterial genomes*)

Five species / eight genomes, hand-checked against the identifiers in the OSDR tables; tables are in `data/bacteria/`.

| OSD | species | locus tags / ids | genome (`GENOMES` key) |
|---|---|---|---|
| 14, 15 | *P. aeruginosa* PAO1 | Affymetrix probesets `PA####_…_at` | `paer_pao1` (GCF_000006765.1) |
| 554 | *P. aeruginosa* PA14 | `PA14_RS#####` | `paer_pa14` (GCF_000014625.1) |
| 138, 185 (215, 216 raw only) | *B. subtilis* 168 | `BSU_#####` | `bsub_168` (GCF_000009045.1) |
| 145 | *S. aureus* UAMS-1 | `QV15_#####` | `saur_uams1` (GCA_000788115.1) |
| 43 | *E. coli* K-12 MG1655 | Affymetrix probesets → Entrez → `b####` | `ecol_k12` (GCF_000005845.2) |
| 95 | *E. coli* ATCC 4157 (NCTC 86) | `E4V11_RS#####` | `ecol_39282` (GCF_900699165.1) |
| 728 (raw only) | *E. coli* B REL606 | `ECB_RS#####` | `ecol_rel606` (GCF_000017985.1) |
| 31, 39 | *R. rubrum* ATCC 11170 | array reporters `RruA####` | `rrub_11170` (GCF_000013085.1) |

```
gene-lingo bacteria build -o data/bacteria        # downloads genomes, writes tables (~1.5 min)
gene-lingo bacteria map --from OSD-554 --to OSD-138 -d data/bacteria PA14_RS00005 PA14_00010
```

`hub.tsv` has one row per *B. subtilis* CDS with an ortholog column per genome, `core.tsv` the 488 genes with an ortholog in all five species, `pairs/<a>__<b>.tsv` the direct hits for all 28 genome pairs, `genes.tsv.gz` every gene with its alternative names, `genomes.json` the genome registry.

### 2. Any genomes, not only OSDR (`gene-lingo genomes`, page *Any genomes*)

Give accessions, organism/strain names, taxids or your own annotation files; any locus-tag convention works because genes are matched by protein sequence.

```
gene-lingo genomes find "Aliivibrio fischeri ES114"        # look up an assembly (RefSeq index, cached after the first download)
gene-lingo genomes build "Bacillus subtilis 168" GCF_000011805.1 taxid:1280 mygenome.gff.gz,mygenome.faa.gz,MyStrain \
    -o out/ [--hub "Bacillus subtilis 168"] [--hub-only]
gene-lingo bacteria map -d out/ --from "B. subtilis 168" --to GCF_000011805.1 BSU_00010
```

Finished ortholog comparisons are cached by assembly accession, so adding a genome only computes the new pairs; `--hub-only` compares each genome with the hub only (n-1 instead of n(n-1)/2) and other pairs are then derived through the hub.

### 3. Any OSDR study (`gene-lingo osdr`, page *Any OSDR study*)

```
gene-lingo osdr info OSD-305 OSD-356 OSD-554       # organism, strain, gene table, and the genome chosen with the evidence
gene-lingo osdr build OSD-305 OSD-356 OSD-554 -o out/
gene-lingo osdr survey -o survey.tsv               # every bacterial transcription study in OSDR
gene-lingo bacteria map -d out/ --from OSD-305 --to OSD-554 M495_RS00005
```

For each study the organism, NCBI taxon and strain come from the OSDR metadata. Candidate genomes are ranked: an assembly named in the metadata, then the strain by name (including a parent strain named in the text, e.g. "derivative of SL1344"), then the species' reference genome. If the study has a processed gene table (counts or microarray; probe ids and Entrez ids are understood) the candidates are checked against its identifiers and the answer is graded *validated*; studies with raw reads only get the best name match, graded *stated* / *strain* / *species*. Override a guess on the page or with `genomes build`.

`data/osdr_bacterial_survey.tsv` is the result for all 28 bacterial transcription studies in OSDR (9 validated against identifiers, 1 stated in the metadata, 14 by strain name, 4 species-level), and `data/osdr_all/` the ortholog tables for the 17 genomes involved (hub-only, *B. subtilis* 168 as hub; 230 hub genes have an ortholog in all 13 species):

| Study | Organism | Gene table | Genome chosen | How sure |
|---|---|---|---|---|
| OSD-11 | *Salmonella enterica* | none processed | S. enterica LT2 (`GCF_000006945.2`) | species-level guess |
| OSD-14 | *Pseudomonas aeruginosa* | 5900 ids | P. aeruginosa PAO1 (`GCF_000006765.1`) | validated (98% of ids found) |
| OSD-15 | *Pseudomonas aeruginosa* | 5900 ids | P. aeruginosa PAO1 (`GCF_000006765.1`) | validated (98% of ids found) |
| OSD-28 | *Bacillus subtilis* | none processed | B. subtilis 168 (`GCF_000009045.1`) | strain name matched |
| OSD-31 | *Rhodospirillum rubrum* | none processed | R. rubrum ATCC 11170 (`GCF_000013085.1`) | species-level guess |
| OSD-43 | *Escherichia coli* | 10208 ids | E. coli K-12 substr. MG1655 (`GCF_000005845.2`) | validated (68% of ids found) |
| OSD-81 | *Bacillus subtilis* | none processed | B. subtilis 168 (`GCF_000009045.1`) | strain name matched |
| OSD-90 | *Mycobacterium marinum* | none processed | M. marinum 1218R (`GCF_003391415.1`) | strain name matched |
| OSD-95 | *Escherichia coli* | 4943 ids | E. coli NCTC86 (`GCF_900699165.1`) | validated (92% of ids found) |
| OSD-106 | *Streptococcus mutans* | none processed | S. mutans UA159 (`GCF_000007465.2`) | strain name matched |
| OSD-119 | *Vibrio fischeri ES114* | none processed | A. fischeri ES114 (`GCF_000011805.1`) | strain name matched |
| OSD-138 | *Bacillus subtilis* | 4443 ids | B. subtilis 168 (`GCF_000009045.1`) | validated (98% of ids found) |
| OSD-145 | *Staphylococcus aureus* | 2808 ids | S. aureus UAMS-1 (`GCA_000788115.1`) | validated (94% of ids found) |
| OSD-185 | *Bacillus subtilis* | 4443 ids | B. subtilis 168 (`GCF_000009045.1`) | validated (98% of ids found) |
| OSD-199 | *Lactobacillus acidophilus* | none processed | L. acidophilus ATCC 4356 (`GCF_034298135.1`) | strain name matched |
| OSD-215 | *Bacillus subtilis* | none processed | B. subtilis 168 (`GCF_000009045.1`) | strain name matched |
| OSD-216 | *Bacillus subtilis* | none processed | B. subtilis 168 (`GCF_000009045.1`) | strain name matched |
| OSD-277 | *Salmonella enterica subsp. enterica serovar Typhimurium* | none processed | S. enterica SL1344 (`GCF_000210855.2`) | strain name matched |
| OSD-305 | *Serratia liquefaciens* | 5006 ids | S. liquefaciens ATCC 27592 (`GCF_000422085.1`) | validated (96% of ids found) |
| OSD-323 | *Salmonella enterica subsp. enterica serovar Typhimurium* | none processed | S. enterica LT2 (`GCF_000006945.2`) | species-level guess |
| OSD-356 | *Vibrio fischeri* | none processed | A. fischeri ES114 (`GCF_000011805.1`) | strain name matched |
| OSD-452 | *Escherichia coli* | none processed | E. coli K-12 substr. MG1655 (`GCF_000005845.2`) | strain name matched |
| OSD-528 | *Mycobacterium marinum* | none processed | M. marinum 1218R (`GCF_003391415.1`) | strain name matched |
| OSD-554 | *Pseudomonas aeruginosa* | 6041 ids | P. aeruginosa UCBPP-PA14 (`GCF_000014625.1`) | validated (98% of ids found) |
| OSD-728 | *Escherichia coli* | none processed | E. coli REL606 (`GCF_000017985.1`) | strain name matched |
| OSD-729 | *Burkholderia contaminans* | none processed | B. contaminans B17-01563-1 (`GCF_022533485.1`) | assembly named in metadata |
| OSD-810 | *Paraburkholderia phymatum* | none processed | P. phymatum STM815 (`GCF_000020045.1`) | species-level guess |
| OSD-987 | *Streptococcus mutans* | none processed | S. mutans UA159 (`GCF_000007465.2`) | strain name matched |

Studies whose genome is graded below *validated* have no processed gene table on OSDR, so the choice cannot be checked; a wrong strain only matters for strain-specific genes.

### Checks and limits

Run on real data (offline unit tests in `tests/`): every protein-coding locus tag in the OSD-95/138/145/185/554 tables resolves to a gene (the rest are tRNA/rRNA/pseudogenes) and 96 % of the OSD-305 identifiers do; PAO1↔PA14 and the three *E. coli* strains give 3.8–5.3k reciprocal hits at ≥99.6 % median identity; for non-hub pairs the direct hit and the route through *B. subtilis* agree for 94.7 % of genes; 96.6 % of the 3,786 *R. rubrum* array oligos lie inside the gene the `RruA####` → `Rru_A####` rule assigns; the general pipeline reproduces the curated *B. subtilis*–PA14 pairs exactly (1153/1153). Reciprocal best hit gives one ortholog per gene, so paralog families can resolve to the wrong member (e.g. *B. subtilis* `lysS` → *E. coli* `lysU`); use the `pident` column to triage. Multi-strain arrays (OSD-43 covers K-12, O157:H7 and CFT073) recognise several genomes, so the strain named in the study wins within a margin.
