import argparse
import numpy as np

from brute_force import brute_force_enhanced
from triangle import trianglex


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--file_path",
        type=str,
        required=True,
        help="file path to the pdb file of the protein of interest",
    )
    parser.add_argument(
        "-a",
        "--algorithm",
        action="store_false",
        help="set flag to use triangle search instead of brute force",
    )
    parser.add_argument(
        "-m",
        "--motif_types",
        type=str,
        required=True,
        help="3letter codes of residues in the motif as string like SER-HIS-ASP",
    )
    parser.add_argument(
        "-d",
        "--dists",
        required=True,
        type=str,
        help="distances between all residues as string like 38.44-26.00-33.31 "
        "distances are like from itertools import combinations "
        "list(combinations(range(len(motif_res)), 2)) which results in "
        "[(0, 1), (0, 2), (1, 2)] what means the first distance is between SER-HIS, "
        "second between SER-ASP and the third between HIS-ASP",
    )
    parser.add_argument(
        "-c",
        "--cushion",
        type=str,
        required=False,
        default=2.0,
        help="max deviation of the given distances between residues as string like "
        "2.0-3.0-4.0 or only one number like 2.0 so the same cushion will be used for "
        "all distances",
    )
    parser.add_argument(
        "-fb",
        "--feedback",
        action="store_false",
        help="set flag do suppress terminal output",
    )
    parser.add_argument(
        "-y",
        "--pymol_check",
        action="store_true",
        help="set flag do get pymol commands for the found motifs as pseudoatoms",
    )
    parser.add_argument(
        "-p",
        "--pdb",
        type=str,
        required=False,
        default=None,
        help="pdb code for data of protein of interest from the pdb website without "
        "saving the file and target_pdb_file will be ignored",
    )

    args = parser.parse_args()
    # convert string inputs to numpy arrays with needed dtype and size
    motif_types_conv = np.asarray(args.motif_types.split("-"))
    dists_conv = np.asarray(args.dists.split("-"), dtype=float)
    num_dists = len(dists_conv)
    if type(args.cushion) == float:
        cushion = np.asarray([2.0] * num_dists)
    else:
        cushion = np.asarray(args.cushion.split("-"), dtype=float)
        if num_dists > 1 and len(cushion) == 1:
            cushion = np.asarray([float(cushion)] * num_dists)

    d = {
        "path_to_prot": args.file_path,
        "motif_types": motif_types_conv,
        "motif_distances": dists_conv,
        "cushion": cushion,
        "feedback": args.feedback,
        "pymol_check": args.pymol_check,
        "from_pdb": args.pdb,
    }
    
    if args.algorithm:
        brute_force_enhanced(**d)
    else:
        trianglex(**d)



if __name__ == "__main__":
    main()
