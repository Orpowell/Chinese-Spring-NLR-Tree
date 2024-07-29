from nlr_genes import get_NLR_genes
from NLR_re_annotator import itol_nbd_nbarc


def nlr_finder(cds: str, annotation: str, proteins: str, email: str, output: str, cloned_nlrs: str = None, motifs: str = None, nbd: bool = False) -> None:

    meerkat = get_NLR_genes(
        nlr_annotations=annotation,
        genome_proteins=proteins,
        genome_cds=cds,
        out=output,
        email=email,
        additional_nlrs=cloned_nlrs
    )

    if nbd and motifs is not None:

        itol_nbd_nbarc(NLRS=annotation, 
                    NLR_MOTIFS=motifs, 
                    proteins=meerkat)