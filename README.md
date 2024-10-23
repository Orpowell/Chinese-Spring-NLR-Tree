# NLR Finder

NLR Finder identifies and extracts the best (longest) protein sequence from transcripts annotated by NLR-annotator for phylogenetic analysis and visualisation with [iTOL](https://itol.embl.de/). Note this tool has only been tested on the Chinese Spring v2.1 Reference Genome and may not work for another genomes (see Issues). 

This code is based on an analysis performed by [Burkhard Steuernagel](https://github.com/steuernb/wheat_nlr/blob/master/java/src/PreparePhylogenetics.java) for [The NLR-Annotator Tool Enables Annotation of the Intracellular Immune Receptor Repertoire](https://doi.org/10.1104/pp.19.01273).

## Dependencies

    Pandas v2.1.4
    Biopython v1.78
    NLR-annotator v2.1 (required for NLR annotations)
    iTOL v6

## Usage

NLR Finder can be run using the following command:

    python3 src/main.py \
    -c cds.fasta \
    -a annotaiton.txt \
    -p protein.fasta \
    -o output.faa \

## Parameters

The following parameters are required to run NLR Finder:

### Required Parameters
| parameter | description |
|---|---|
| -c, --cds | Coding sequences from a reference genome |
| -a, --annotation | NLR-Annotator output (-o output.txt) of coding sequences |
| -p --protein | Sequences corresponding to coding sequences |
| -o, --output | Name of the output fasta file |


### Optional Parameters

The following optional parameters can be provided to add cloned sequences and annotate NBD-NBARC NLRs: 

| parameter | description |
|---|---|
| -x, --cloned_nlrs | 2 column headerless tsv file of gene name (used for annotaion) and GenBank accession number (used to fetch sequence). An example can be found in example.tsv |
| -e, --email | required for accessing NCBI API (required when using --cloned_nlrs) |
| -nbd | Enables annotation of NBD-NBARC NLRs |
| -m, --motifs| Motifs from NLR-Annotator (-m output.motifs.bed) for CDS (required when using -nbd) |

## Outputs

| output | description |
|---|---|
| output.faa | A fasta file of NLR protein sequence identified from provided cds and sequences specified by --cloned_nlrs. File name specified by --output |
| itol_nlr_labels.txt | iTOL annotations for sequences provided by --cloned_nlrs. The corresponding sequences in tree are labelled with the name provided in the tsv.|
| itol_nbs_nbarc_stars.txt | iTOL annotations for sequences re-annotated as NBD-NBARC NLRs. The corresponding sequences in tree are denoated by a red star. |

## Using NLR Finder to build a Phyolgenetic Tree

The following is a mock analysis using NLR Finder to create a phylogenetic tree of unique NLRs from a reference genome:

1. Run NLR-Annotator on reference cds

2. Run NLR Finder

3. Generate multiple sequence alignment from NLR Finder output

4. Generate Phylogenetic Tree

5. Load tree into [iTOL](https://itol.embl.de) and visualise annotations from NLR Finder


## Questions and Issues

Any issues can be reported on the repositories issues page. Alternatively, you can reach out and contact me directly [here](mailto:nlr-finder@oliverpowell.com).

This tool was designed to analyse the Chinese Spring v2.1 Reference Genome and may not work with other genomes. If you are interested in applying this analysis to your genome of interest please contact me and I'd be happy to help!

## Citing this tool

If you use this tool in your research please cite:

- [Steuernagel et al.: The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire, Plant Physiology, 2020](https://www.ncbi.nlm.nih.gov/pubmed/32184345)

- [Chen et al.: A wheat tandem kinase sensor activates an NLR helper to trigger immunity, 2024](https://doi.org/10.1101/2024.08.30.610287)


