from nlr_genes import get_NLR_genes
from NLR_re_annotator import itol_nbd_nbarc

# Identify NLR genes
CS_CDS_NLR_ANNOTATION = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/analysis/nlr_annotation/CDS/CS_CDS_NLR_annotator.txt"
CS_CDS_NLR_MOTIFS = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/analysis/nlr_annotation/CDS/CS_CDS_NLR_annotator.motifs.bed"
CS_PROTEIN_FASTA = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/data/chinese_spring/protein.faa"
CS_CDS_FASTA = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/data/chinese_spring/cds_from_genomic.fna"
NLR_DOMAIN_ANNOTATION = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/analysis/NLR_sequence_analysis/pfam_annotation.txt"
EMAIL = "oliver.powell@kaust.edu.sa"
added_nlrs = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/data/additional_nlr/selected_nlrs.tsv"

# Get all NLR proteins from CDS (longest one if multiple transcripts)
meerkat = get_NLR_genes(
    nlr_annotations=CS_CDS_NLR_ANNOTATION,
    genome_proteins=CS_PROTEIN_FASTA,
    genome_cds=CS_CDS_FASTA,
    out="cs_and_cloned_nlrs.fasta",
    email=EMAIL,
    additional_nlrs=added_nlrs
)

itol_nbd_nbarc(NLRS=CS_CDS_NLR_ANNOTATION, 
               NLR_MOTIFS=CS_CDS_NLR_MOTIFS, 
               proteins=meerkat)