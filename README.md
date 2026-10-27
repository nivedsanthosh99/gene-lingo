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

The sidebar has a **How it works** page with the full flowchart, API calls and formats.

Pick input/output species and ID types, paste genes or upload a CSV, and download the results table as CSV. Same species = ID conversion; different species = ortholog lookup, with the orthologs shown in your chosen output ID type.

## Bacterial locus-tag mapping (OSDR studies)

The OSDR bacterial studies use locus tags from different annotations (`PA14_RS00005` vs old `PA14_00010` vs PAO1 `PA0001`, `BSU_00010` vs `BSU00010`, `RRU_RS00005` vs `Rru_A0001`, ...), so name-based lookups (and tools keyed to one annotation) miss many genes. `gene_lingo.bacteria` instead:

1. pins each OSD to the genome its locus tags come from (`GENOMES`; NCBI GFF + protein FASTA),
2. computes reciprocal best DIAMOND hits (both proteins ≥50 % covered, e-value ≤1e-5) between every pair of genomes and against the hub, *B. subtilis* 168,
3. resolves old/alternative locus-tag spellings and array probe ids (OSD-14/15 `PA####_at`, OSD-43 Entrez ids, OSD-31/39 `RruA####`) to the canonical tag.

| OSD | species | locus tags / ids | genome (`GENOMES` key) |
|---|---|---|---|
| 14, 15 | *P. aeruginosa* PAO1 | Affymetrix probesets `PA####_…_at` | `paer_pao1` (GCF_000006765.1) |
| 554 | *P. aeruginosa* PA14 | `PA14_RS#####` | `paer_pa14` (GCF_000014625.1) |
| 138, 185 (215, 216 raw only) | *B. subtilis* 168 | `BSU_#####` | `bsub_168` (GCF_000009045.1) |
| 145 | *S. aureus* UAMS-1 | `QV15_#####` | `saur_uams1` (GCA_000788115.1) |
| 43 | *E. coli* K-12 MG1655 | Affymetrix probesets → Entrez → `b####` | `ecol_k12` (GCF_000005845.2) |
| 95 | *E. coli* 39282_E01-4 | `E4V11_RS#####` | `ecol_39282` (GCF_900699165.1) |
| 728 (raw only) | *E. coli* B REL606 | `ECB_RS#####` | `ecol_rel606` (GCF_000017985.1) |
| 31, 39 | *R. rubrum* ATCC 11170 | array reporters `RruA####` | `rrub_11170` (GCF_000013085.1) |

```
sudo apt install diamond-aligner
gene-lingo bacteria build -o data/bacteria            # downloads genomes, writes tables (~1.5 min)
gene-lingo bacteria map --from OSD-554 --to OSD-138 -d data/bacteria PA14_RS00005 PA14_00010
```

Output in `data/bacteria/`: `hub.tsv` (one row per *B. subtilis* CDS, one ortholog column per genome), `core.tsv` (488 *B. subtilis* genes with an ortholog in all five species), `pairs/<a>__<b>.tsv` (direct reciprocal best hits for all 28 genome pairs).

Checks run on the real data (see `tests/test_bacteria.py` for the offline unit tests): every protein-coding locus tag in the OSD-95/138/145/185/554 count tables resolves to a CDS (the rest are tRNA/rRNA/pseudogenes); PAO1↔PA14 and the three *E. coli* strains give 3.8–5.3k reciprocal hits at ≥99.6 % median identity; for non-hub pairs the direct hit and the route through *B. subtilis* agree for 94.7 % of genes; 96.6 % of the 3,786 *R. rubrum* array oligos lie inside the CDS of the locus tag the `RruA####` → `Rru_A####` rule assigns. Known limits: reciprocal best hit picks one ortholog per gene, so paralog families can resolve to the wrong member (e.g. *B. subtilis* `lysS` → *E. coli* `lysU`); use the `pident` column in `pairs/` to triage. OSD-215/216/728 have no processed gene tables yet; their genome choice is inferred from the study metadata.
