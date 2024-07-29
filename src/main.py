import os
import sys
import logging
import argparse
from nlr_finder import nlr_finder


logging.basicConfig(
    stream=sys.stderr,
    format="%(asctime)s - %(message)s",
    datefmt="%d-%b-%y %H:%M:%S",
    level=logging.INFO,
)


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error(f"Input file ({arg}) not found!")

    else:
        return arg

def main():
    # Initialise parser and subparser for each command
    parser = argparse.ArgumentParser()

    # Blast2bed arguments
    parser.add_argument(
        "-c",
        "--cds",
        metavar="cds",
        required=True,
        type=lambda x: is_valid_file(parser, x),
        help="CDS of interest",
    )
    
    parser.add_argument(
        "-a",
        "--annotation",
        metavar="annotation",
        required=True,
        type=lambda x: is_valid_file(parser, x),
        help="NLR-Annotator annotations for CDS",
    )

    parser.add_argument(
        "-p",
        "--proteins",
        metavar="proteins",
        required=True,
        type=lambda x: is_valid_file(parser, x),
        help="Protein sequences corresponding to the CDS",
    )

    parser.add_argument(
        "-e",
        "--email",
        metavar="email",
        required=True,
        type=str,
        help="Email. Required for access to NCBI",
    )

    parser.add_argument(
        "-o",
        "--output",
        metavar="output",
        required=True,
        type=str,
        help="Output file name",
    )

    parser.add_argument(
        "-x",
        "--cloned_nlrs",
        metavar="cloned_nlrs",
        required=False,
        default=None,
        type=lambda x: is_valid_file(parser, x),
        help="TSV with gene label and genbank id",
    )

    parser.add_argument(
        "-m",
        "--motifs",
        metavar="motifs",
        required=False,
        default=None,
        type=lambda x: is_valid_file(parser, x),
        help="Motifs from NLR-Annotator",
    )
    parser.add_argument(
        "-nbd",
        required=False,
        default=False,
        action="store_true",
        help="Identify and Annotate NBD-NBARC NLRs",
    )

    args = parser.parse_args()

    nlr_finder(cds=args.cds,
           annotation=args.annotation,
           proteins=args.proteins,
           email=args.email,
           output=args.output,
           cloned_nlrs=args.cloned_nlrs,
           motifs=args.motifs,
           nbd=args.nbd)

if __name__ == "__main__":
    main()