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
