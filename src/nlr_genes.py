import pandas as pd
from Bio import SeqIO
from Bio import Entrez
from Bio import GenBank

pd.options.mode.copy_on_write = True


def get_NLR_transcripts_simple(nlr_annotation: str) -> set[str]:
    """
    :param nlr_annotation: NLR annotations from CDS of genome
    :return
    """

    nlr_df = pd.read_csv(
        nlr_annotation,
        sep="\t",
        header=None,
        names=["cds", "nlr", "domains", "start", "end", "strand", "motifs"],
    )

    return set(nlr_df["cds"].to_list())


def best_transcript_per_gene(cds: str, nlr_transcripts: set[str]) -> set[str]:
    def get_protein_header(header: str) -> str:
        row = header.split("_")

        protein = f"{row[3]}_{row[4]}"

        return protein

    gene_transcript = dict()
    gene_size = dict()

    for record in SeqIO.parse(cds, "fasta"):
        if record.id in nlr_transcripts:
            gene = record.description.split(" ")[1]

            if (gene not in gene_transcript) or (gene_size[gene] < len(record.seq)):
                gene_transcript[gene] = record.id
                gene_size[gene] = len(record.seq)

    return set([get_protein_header(gene) for gene in gene_transcript.values()])


def find_nlr_protein_seqs(
    proteins: str, nlr_transcripts: set[str], output: str
) -> None:
    """
    :param proteins: path to protein sequences for NLRs
    :param nlr_transcripts: set of all NLR transcript IDs
    :return None: NLR protein sequences are written to a file
    """

    nlrs = []

    nlrs_added = 0
    for record in SeqIO.parse(proteins, "fasta"):
        if record.id in nlr_transcripts:
            nlrs.append(record)
            nlrs_added += 1

    SeqIO.write(nlrs, output, "fasta")

    print(f"{nlrs_added=}")


def get_NLR_genes(
    nlr_annotations: str, genome_cds: str, genome_proteins: str, out: str, email: str, additional_nlrs=None 
) -> set[str]:
    # nlr_transcripts = get_NLR_transcripts(CS_CDS_NLR_ANNOTATION)
    meerkat_transcripts = get_NLR_transcripts_simple(nlr_annotation=nlr_annotations)

    # all_nlr_genes = best_transcript_per_gene(CS_CDS_FASTA, nlr_transcripts)
    meerkat_genes = best_transcript_per_gene(genome_cds, meerkat_transcripts)

    # find_nlr_protein_seqs(CS_PROTEIN_FASTA, all_nlr_genes, "NLR_proteins.best_gene.fasta")
    find_nlr_protein_seqs(genome_proteins, meerkat_genes, out)

    if additional_nlrs is not None:
        cloned: list[cloned_NLR] = fetch_additional_nlrs(added_nlrs = additional_nlrs, email=email)
        write_itol_labels(cloned)
        write_nlrs(nlrs = cloned, out=out)

    return meerkat_genes


class cloned_NLR:
    def __init__(self, name, genbank, sequence) -> None:
        self.name = name
        self.genbank = genbank
        self.sequence = sequence

    def format_fasta(self) -> str:
        return f">{self.genbank}\n{self.sequence}\n"
    
    def format_itol_label(self) -> str:
        return f"{self.genbank},{self.name},-1,#000000,italic,10,0\n"


def fetch_additional_nlrs(added_nlrs: str, email: str) -> list[cloned_NLR]:
    Entrez.email = email

    df = pd.read_csv(added_nlrs, sep="\t", header=None, index_col=None)

    nlr_genebank = dict(zip(df[0], df[1]))
    genebank_nlr = dict(zip(df[1], df[0]))

    genbank_ids = ",".join(list(nlr_genebank.values()))

    handle = Entrez.efetch(db="protein", id=genbank_ids, rettype="gb", retmode="text")

    records = GenBank.parse(handle)

    cloned_nlrs = [cloned_NLR(name = genebank_nlr[record.accession[0]], genbank=record.accession[0],sequence=record.sequence) for record in records]
        
    return cloned_nlrs


def write_itol_labels(nlrs: list[cloned_NLR]) -> None:
    with open("itol_nlr_labels.txt", "w+") as file:
        file.write('DATASET_TEXT\n')
        file.write('SEPARATOR COMMA\n')
        file.write('DATASET_LABEL,NLR labels\n')
        file.write("COLOR,#ff0000\n")
        file.write("MARGIN,0\n\n\n")
        file.write("DATA\n")

        for nlr in nlrs:
            file.write(nlr.format_itol_label())


def write_nlrs(nlrs: list[cloned_NLR], out) -> None:
    with open(out, "a+") as file:
        for nlr in nlrs:
            file.write(nlr.format_fasta())
