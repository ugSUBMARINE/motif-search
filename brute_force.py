import itertools
from itertools import combinations

import numpy as np

from ms_utils import pseudoatom_positions, X_DICT, pymol_print


def brute_force_enhanced(
    path_to_prot: str,
    motif_types: list[str],
    motif_distances: np.ndarray[tuple[int, int], np.dtype[int | float]],
    cushion: np.ndarray[tuple[int], np.dtype[int | float]],
    feedback: bool = False,
    pymol_check=False,
    from_pdb: str | None = None,
):
    """*** This function was writen by Karl Gruber ***
    finds active site motifs in a given pdb file
    :parameter
        - path_to_prot:
          file path to pdb file
        - motif_types:
          3letter codes of residues in the motif like [SER, HIS, ASP]
        - motif_distances:
          distances between all residues eg [38.44 26.00 33.31]
          distances are like
          from itertools import combinations
          list(combinations(range(len(motif_types)), 2)) which results in
          [(0, 1), (0, 2), (1, 2)]
          what means the first distance is between SER-HIS, second between SER-ASP and
          the third between HIS-ASP
        - cushion:
          max deviation of the given distances between residues
        - feedback
          if intermediate output should be written in terminal
        - pymol_check:
          if True the coordinates of the pseudoatom positions of the valid residues
          get printed in a way that one can copy and past them in the pymol terminal
          and get them displayed and the valid residues get shown as sticks
        - from_pdb:
          pdb code for data of protein of interest from the pdb website
          without saving the file and target_pdb_file will be ignored

    :return
        - res_data:
          each motif is in a 2D array like
          [[3letter code, chain id, residue number],...]
        - res_coords:
          xyz coordinates of each residue of a motif is in a 2D array like
          [[x1,y1,z1],...]
    """
    n_motif = len(motif_types)
    n_dist = len(motif_distances)
    needed_dists = n_motif * (n_motif - 1) / 2
    print(needed_dists)
    if n_dist != needed_dists:
        raise ValueError(
            f"Wrong number of distances {n_dist} instead of {needed_dists}"
        )
    else:
        if feedback:
            print("EVERYTHING FINE!")

    res_data, res_coords = pseudoatom_positions(
        atom_pos_dict=X_DICT, target_pdb_file=path_to_prot, pdb_code=from_pdb
    )
    res_data = np.array(res_data)
    resn = res_data[:, 0]
    chain = list(res_data[:, 1])
    resi = list(res_data[:, 2])

    # calculate all pairwise distances
    coords = np.array(res_coords)
    dist = np.sqrt(
        np.sum((coords.reshape(-1, 1, 3) - coords.reshape(1, -1, 3)) ** 2, axis=2)
    )

    # extract indices of residues that are potentially part in the
    n_motifs = 1
    res_indices = []
    for r in motif_types:
        indices = np.where(resn == r)[0]  # indices of all residues of a certain type
        n_motifs *= len(indices)
        res_indices.append(indices)
    if feedback:
        print("Number of possible combinations: {}".format(n_motifs))
        for i, indices in enumerate(res_indices):
            print("{}: {}".format(motif_types[i], indices))
        print()

    def analyze_pairs(res_ind_1, res_ind_2, distances, ideal_dist, cushion):
        """
        returns lists of indices for residue pairs which obey the distance criterion
        """
        res_1 = []
        res_2 = []
        for pair in itertools.product(res_ind_1, res_ind_2):
            if pair[0] != pair[1]:  # must be different residues
                if abs(distances[pair] - ideal_dist) <= cushion:
                    res_1.append(pair[0])
                    res_2.append(pair[1])

        return res_1, res_2

    n_3 = 0
    for n, (i, j) in enumerate(itertools.combinations(range(n_motif), 2)):
        if feedback:
            print("looking at {}-{} pairs".format(motif_types[i], motif_types[j]))
            print(
                "distance should be {:.2f} +/- {:.2f}".format(
                    motif_distances[n], cushion[n]
                )
            )

        # get those residues that obey the distance criterion
        n_3 += len(res_indices[i]) * len(res_indices[j])
        res_1, res_2 = analyze_pairs(
            res_indices[i], res_indices[j], dist, motif_distances[n], cushion[n]
        )
        if feedback:
            print("{} pairs found".format(len(res_1)))

        res_1 = list(set(res_1))
        res_2 = list(set(res_2))
        if feedback:
            print("remaining {} residues: {}".format(motif_types[i], res_1))
            print("remaining {} residues: {}".format(motif_types[j], res_2))
            print()

        # modify the residue list for the next iteration of the loop
        res_indices[i] = np.array(res_1)
        res_indices[j] = np.array(res_2)

    if feedback:
        print("{} pairs checked\n".format(n_3))
    if feedback:
        print("reduced res_indices")
        print(res_indices)
    n_motifs_reduced = 1
    for ind in res_indices:
        n_motifs_reduced *= len(ind)
    if feedback:
        print(n_motifs_reduced, " vs. ", n_motifs)
        print()

    n_1 = 0
    n_2 = 0
    all_found_motives = []
    for motif in itertools.product(*res_indices):
        if len(motif) == len(set(motif)):  # use only combinations with unique indices
            n_1 += 1
            dist_ind = tuple(zip(*itertools.combinations(motif, 2)))
            d = dist[dist_ind]
            flag = np.all(np.abs(d - motif_distances) <= cushion)
            if flag:
                all_found_motives += [list(motif)]
                if feedback:
                    print("Motif found:")
                    for i, res in enumerate(motif):
                        print(
                            "{:4d}: {:3s}-{:1s}{:4s}".format(
                                i + 1, resn[res], chain[res], resi[res]
                            )
                        )
                    print("Distances:")
                    print(d)
                    print("Absolute distances from ideality:")
                    print(np.abs(d - motif_distances))
                    print()
        else:
            n_2 += 1
    if feedback:
        print("Unique motifs: {}".format(n_1))
        print("Discarded non-unique motifs: {}".format(n_2))

    result_indices = np.unique(np.sort(all_found_motives), axis=0)
    if len(result_indices) > 0:
        res_data = np.asarray(res_data)
        res_coords = np.asarray(res_coords)
        if pymol_check:
            pymol_print(res_data, res_coords)
        return (
            res_data[result_indices],
            res_coords[result_indices],
        )
    else:
        if feedback:
            print("+++ Motif not found +++")
        return


if __name__ == "__main__":
    brute_force_enhanced(
        path_to_prot="/PATH/TO/PDB/FILE",
        motif_types=["SER", "ASN", "ARG", "SER"],
        motif_distances=np.asarray([7.31146104, 7.84910799, 7.72271787, 9.0, 9.0, 9.0]),
        cushion=np.asarray([2] * 6),
        feedback=True,
        pymol_check=True,
    )
