from nlr_finder import nlr_finder

# Identify NLR genes
CS_CDS_NLR_ANNOTATION = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/analysis/nlr_annotation/CDS/CS_CDS_NLR_annotator.txt"
CS_CDS_NLR_MOTIFS = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/analysis/nlr_annotation/CDS/CS_CDS_NLR_annotator.motifs.bed"
CS_PROTEIN_FASTA = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/data/chinese_spring/protein.faa"
CS_CDS_FASTA = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/data/chinese_spring/cds_from_genomic.fna"
EMAIL = "oliver.powell@kaust.edu.sa"
OUTPUT = "CS_and_cloned.test"
added_nlrs = "/home/powellor/Documents/projects/sr62_renjie/nlr_phylogenetics/data/additional_nlr/selected_nlrs.tsv"

# Get all NLR proteins from CDS (longest one if multiple transcripts)
nlr_finder(cds=CS_CDS_FASTA,
           annotation=CS_CDS_NLR_ANNOTATION,
           proteins=CS_PROTEIN_FASTA,
           email=EMAIL,
           output=OUTPUT,
           cloned_nlrs=added_nlrs,
           motifs=CS_CDS_NLR_MOTIFS,
           nbd=True)