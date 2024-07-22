import pandas as pd

pd.options.mode.copy_on_write = True


class NLR:
    def __init__(self, name, chromosome, start, end, strand, annotation) -> None:
        self.name: str = name
        self.chromosome: str = chromosome
        self.start: int = start
        self.end: int = end
        self.strand: str = strand
        self.annotation: str = annotation
        self.nbs_nbarc: bool = False
        self.gene_loci: list[str] = []
        self.protein: str = "_".join(name.split("_")[3:5])

    def get_motifs(self, motif_df: pd.DataFrame) -> None:
        valid_motifs = motif_df[
            (motif_df[0] == self.chromosome)
            & (motif_df[1] >= self.start)
            & (motif_df[2] <= self.end)
            & (motif_df[5] == self.strand)
        ]

        self.motifs = valid_motifs[3].to_list()

    def check_nbs_nbarc(self) -> None:
        m1: int = self.motifs.count("motif_1")
        m4: int = self.motifs.count("motif_4")
        m5: int = self.motifs.count("motif_5")
        m6: int = self.motifs.count("motif_6")

        nbarc_motif_counts: list[int] = [m1, m4, m5, m6]

        if list(map(lambda x: x > 1, nbarc_motif_counts)).count(True) > 1:
            self.nbs_nbarc = True

    def get_gene_loci(self, gff_df: pd.DataFrame) -> None:
        valid_loci = gff_df[
            (gff_df[0] == self.chromosome)
            & (
                ((gff_df[3] <= self.end) & (self.start <= gff_df[4]))
                | ((self.end <= gff_df[3]) & (gff_df[4] <= self.start))
            )
            & (gff_df[6] == self.strand)
        ]

        self.gene_loci = valid_loci.locus.to_list()

    def run(self, motifs) -> None:
        self.get_motifs(motifs)
        self.check_nbs_nbarc()


def write_itol_star_nbd_nbarcs(nbs_nbarc_nlrs):
    with open("itol_nbs_nbarc_stars.txt", "w+") as output:
        output.write("DATASET_BINARY\n")
        output.write("SEPARATOR COMMA\n")
        output.write("DATASET_LABEL,Kinase Domains\n")
        output.write("COLOR,#00ff00\n")
        output.write("FIELD_SHAPES,3\n")
        output.write("FIELD_LABELS,kinase_present\n")
        output.write("FIELD_COLORS,#ff0000\n")
        output.write("DATA\n")

        for protein in nbs_nbarc_nlrs:
            identifiers = protein.name.split("_")
            name = f"{identifiers[3]}_{identifiers[4]}"
            output.write(f"{name},1\n")

def write_itol_range_nbd_nbarcs(nbs_nbarc_nlrs):
    with open("itol_nbs_nbarc_range.txt", "w+") as output:
        output.write("DATASET_RANGE\n")
        output.write("SEPARATOR COMMA\n")
        output.write("DATASET_LABEL,Kinase Domains\n")
        output.write("COLOR,#00ff00\n")
        output.write("FIELD_SHAPES,3\n")
        output.write("FIELD_LABELS,kinase_present\n")
        output.write("FIELD_COLORS,#ff0000\n")
        output.write("DATA\n")

        names = [n.protein for n in nbs_nbarc_nlrs]
            
        clade: str = "|".join(names)
        output.write(f"{clade},{clade},#FF0000\n")

def itol_nbd_nbarc(NLRS: str, NLR_MOTIFS: str, proteins: set[str]):
    all_motifs = pd.read_csv(NLR_MOTIFS, skiprows=2, header=None, sep="\t")

    NLR_array = []

    with open(NLRS) as nlrs:
        for line in nlrs:
            row = line.split("\t")

            NLR_array.append(
                NLR(
                    name=row[1],
                    chromosome=row[0],
                    start=int(row[3]),
                    end=int(row[4]),
                    strand=row[5],
                    annotation=row[2],
                )
            )

    unique_nlrs = [nlr for nlr in NLR_array if nlr.protein in proteins]
        
    [n.run(all_motifs) for n in unique_nlrs]

    nbs_nbarc_nlrs = [nlr for nlr in unique_nlrs if nlr.nbs_nbarc]

    nbs_nbarc_count = len(nbs_nbarc_nlrs)

    print(f"NLRs Identified = {nbs_nbarc_count}")

    write_itol_star_nbd_nbarcs(nbs_nbarc_nlrs=nbs_nbarc_nlrs)