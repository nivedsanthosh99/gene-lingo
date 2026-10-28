import gzip
import random

from gene_lingo import bacteria

AA = "ACDEFGHIKLMNPQRSTVWY"


def prot_seq(rng, n=150):
    return "".join(rng.choice(AA) for _ in range(n))


def mutate(rng, s, frac):
    return "".join(rng.choice(AA) if rng.random() < frac else c for c in s)


def make_genome(tmp_path, key, genes):
    """genes: [(locus_tag, old_tag, symbol, seq)] -> Proteome built from a synthetic GFF + FASTA."""
    gff, faa = tmp_path / f"{key}.gff.gz", tmp_path / f"{key}.faa.gz"
    with gzip.open(gff, "wt") as g, gzip.open(faa, "wt") as f:
        g.write("##gff-version 3\n")
        for i, (tag, old, sym, seq) in enumerate(genes):
            attrs = f"ID=gene-{tag};Name={sym};gene={sym};gbkey=Gene;locus_tag={tag}"
            if old:
                attrs += f";old_locus_tag={old}"
            attrs += f";Dbxref=GeneID:{1000 + i}"
            g.write(f"chr\tRefSeq\tgene\t{i*1000+1}\t{i*1000+450}\t.\t+\t.\t{attrs}\n")
            g.write(f"chr\tRefSeq\tCDS\t{i*1000+1}\t{i*1000+450}\t.\t+\t0\tID=cds-{tag};Parent=gene-{tag};"
                    f"locus_tag={tag};protein_id=P_{tag};product=prod%20{sym}\n")
            f.write(f">P_{tag} prod {sym}\n{seq}\n")
        g.write(f"chr\tRefSeq\tpseudogene\t90000\t90300\t.\t+\t.\tID=gene-PSEUDO;locus_tag=PSEUDO\n")
    return bacteria.Proteome(key, gff, faa)
