import sys
import random
import itertools
from itertools import combinations

import numpy as np

from ms_utils import (
    X_DICT,
    read_motif_dict,
    motif_dict_append,
    pymol_print,
    pseudoatom_positions,
    motif_extraction,
)


def trianglex(
    path_to_prot: str,
    motif_types: list[str],
    motif_distances: list[int | float],
    cushion: int
    | float
    | list[int | float]
    | np.ndarray[tuple[int], np.dtype[int | float]] = 0.5,
    feedback=True,
    pymol_check=False,
    from_pdb: str | None = None,
):
    """searches for motifs defined by motif_types and motif_distances by calculating
        all distances of all atoms against each other, then finding all pairs of
        all residues in motif_types and checking whether they are as close or closer than
        defined in motif_distances - with these pairs - all possible triangles get build -
        after that the motif gets build out of a combination of as many triangles as
        needed to get enough points to build the motif
    :parameter
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
          how much bigger the distances can be and still be valid as found motif
        - feedback:
          if True no internal error messages get printed
        - pymol_check:
          if True the coordinates of the pseudoatom positions of the valid residues
          get printed in a way that one can copy and past them in the pymol terminal
          and get them displayed and the valid residues get shown as sticks
        - from_pdb:
          pdb code for data of protein of interest from the pdb website
          without saving the file and target_pdb_file will be ignored
    :return
        - res:
          each motif is in a 2D array like
          [[3letter code, chain id, residue number],...]
        - coords:
          xyz coordinates of each residue of a motif is in a 2D array like
          [[x1,y1,z1],...]
    """
    # Error handling
    # checking for matching residue and distance inputs
    num_residues = len(motif_types)
    num_dists = len(motif_distances)
    if num_dists != (int(num_residues * (num_residues - 1) / 2)):
        raise ValueError(
            "motif residue number and number of distances don't "
            "define the motif correctly\n"
            f"Number of residues given: {num_residues}\n"
            f"Expeted distances: {int((num_residues * (num_residues - 1)) / 2)}\n"
            f"Given distances: {num_dists}\n"
        )

    # checking cushion and setting cushion_mode for pairwise distance search
    if type(cushion) == float or type(cushion) == int:
        cushion_mode = 0
    elif type(cushion) == list or type(cushion) == np.ndarray:
        if len(cushion) == num_dists:
            cushion_mode = 1
        else:
            raise ValueError(
                "number of distances doesn't match the number of values given "
                "in cushion"
            )
    else:
        raise TypeError(
            f"unsupported operand type(s) - got {type(cushion)} but "
            "'float', 'int' or 'list' expected for cushion"
        )

    # [[3letter, chain id, residue id],...] and [[x1,y1,z1],...] for ever amino acid
    data_and_coords = pseudoatom_positions(X_DICT, path_to_prot, from_pdb)
    res_data = np.asarray(data_and_coords[0])
    res_coords = np.asarray(data_and_coords[1], dtype=float)

    # reshaping for the distance calculation
    coords_a0 = res_coords.reshape(res_coords.shape[0], 1, 3)
    coords_a1 = res_coords.reshape(1, res_coords.shape[0], 3)
    # distances between each point
    distance = np.sqrt(np.sum((coords_a0 - coords_a1) ** 2, axis=2))
    # all possible distance combinations
    combs = list(combinations(range(len(motif_types)), 2))

    # list of lists to store the comb results - indices of residue forming pairs
    valid_pairs = []
    for dist_count, i in enumerate(combs):
        # where aa from motif_types i[0] and i[1] from combs is closer that dist_th
        if cushion_mode == 0:
            cushion_to_use = cushion
        else:
            cushion_to_use = cushion[dist_count]
        # distances allowed are motif_distances[dist_count] +/-cushion
        dist_check = np.where(
            (distance > motif_distances[dist_count] - cushion_to_use)
            & (distance < motif_distances[dist_count] + cushion_to_use)
        )
        val_res1 = dist_check[0]
        val_res2 = dist_check[1]

        # where motif_types i are located in val_res1 and val_res2 and
        # where they form a pair
        pair_test = (res_data[val_res1][:, 0] == motif_types[i[0]]) & (
            res_data[val_res2][:, 0] == motif_types[i[1]]
        )

        # valid residue pairs as their ture index in res_data
        val_pairs = np.stack((val_res1[pair_test], val_res2[pair_test]), axis=1)
        # whether distance comparisons found some matching residues
        if len(val_pairs) == 0:
            if feedback:
                print(
                    f"Error - no {motif_types[i[0]]} + {motif_types[i[1]]} pair found"
                )
            return
        valid_pairs.append(val_pairs.tolist())

    # return for dyads
    if len(motif_types) == 2:
        val_data = res_data[valid_pairs[0]]
        val_coords = res_coords[valid_pairs[0]]
        if pymol_check:
            pymol_print(val_data, val_coords)
        return val_data, val_coords

    def create_check_arr(unequal_arr):
        """combine all pair arrays into one to eliminate duplicates"""
        true_arr = []
        for f in unequal_arr:
            true_arr += f
        return np.unique(np.sort(true_arr), axis=0)

    def create_unique_combinations(
        base_arr,
        check_arr: np.ndarray[tuple[int, int], np.dtype[float]],
        size_penalty: int,
        bool_comb_size: int,
        tri_needed: int | None = None,
    ) -> np.ndarray[tuple[int, int], np.dtype[float]]:
        """creates all possible unique combinations of entries in base_arr,
           where check_arr is used to check whether all combinations of the newly
           build shape feature the subshapes that were previously confirmed - if so
           the new shape is a valid new shape
        :parameter
            - base_arr:
              points to create the new shape
            - check_arr:
              array of subshapes
            - size penalty:
              size (points) of the shape of interest
            - bool_comp_size:
              size of the previous confirmed shape to get all combinations
              from the points in the new shape which then can be checked whether all
              of them were previously confirmed
            - tri_needed:
              number of triangles needed to have enough points to build the motif
        :return
            array with sorted new shapes in a 2D array
        """
        # product for pairs and combinations for array that features unique triangles
        if bool_comb_size == 2:
            to_test = list(itertools.product(*base_arr))
        else:
            to_test = list(itertools.combinations(base_arr, tri_needed))
        shape_storage = []
        # create all possible shapes formed by base arr
        for k in to_test:
            # flatten the combination and produce
            # a list with only the unique residues in the shape k
            shape = list(set(list(itertools.chain.from_iterable(k))))
            # so only the shapes that are the shapes of interest survive
            if len(shape) == size_penalty:
                # whether all combinations that are formed by that shape
                # are really present(check_arr)/possible
                # so all sides of the created shape are somewhere found
                # in valid_pairs
                if shorten_loop(shape, bool_comb_size, check_arr):
                    shape_storage += [shape]
        return np.unique(np.sort(shape_storage), axis=0)

    def shorten_loop(
        sl_shape: int,
        sl_bcs: int,
        sl_check_arr: np.ndarray[tuple[int, int], np.dtype[float]],
    ) -> bool:
        """shortens the loop so when one combination is False it returns False
        and doesn't unnecessarily check the remaining combinations
        :parameter
            - sl_shape:
              shape from create_unique_combinations
            - sl_bsc:
              bool_comb_size from create_unique_combinations
            - sl_check_arr:
              check_arr from create_unique_combinations
        :return
            - bool
        """
        for n in list(itertools.combinations(sl_shape, sl_bcs)):
            if not np.any(np.all(sl_check_arr == np.sort(list(n)), axis=1)):
                return False
        return True

    # all valid formed triangles
    all_tri = []
    # all possible pairs as 2D list to check whether tri contains only matching pairs
    check_pairs = create_check_arr(valid_pairs)
    # combinations of residues that form a triangle -
    # all triangles that can be formed with the residues
    for i in list(combinations(range(len(motif_types)), 3)):
        # pairs of points that are formed inside a triangle
        tri_combs = list(combinations(i, 2))
        pairs_of_tri = []
        for j in tri_combs:
            # where the pairs in valid_pairs can be found
            in_all_combs = np.where(np.all(np.isin(combs, j), axis=1))[0][0]
            tri = np.unique(np.sort(valid_pairs[in_all_combs]), axis=0).tolist()
            pairs_of_tri.append(tri)
        created_tri = create_unique_combinations(pairs_of_tri, check_pairs, 3, 2)
        all_tri += created_tri.tolist()
    # to later check whether the motifs contains all matching triangles
    check_triangle = np.unique(np.sort(all_tri), axis=0)

    # number of triangles required to have enough residues to build the motif
    req_tri = np.ceil(len(motif_types) / 3).astype(int)
    # all created motif
    motif = create_unique_combinations(
        check_triangle, check_triangle, len(motif_types), 3, req_tri
    )

    valid_combinations = np.sort(np.unique(motif, axis=0))
    # checking whether each pair of each valid_combinations
    # is found in a different pair list from valid_pairs
    vc_b = []
    for ci, i in enumerate(valid_combinations):
        incl = []
        # every possible pair combination of residues in i
        for j in np.sort(list(combinations(i, 2))):
            j_sort = np.sort(j)
            # in which valid_pairs pair list j is present
            for cl, l in enumerate(valid_pairs):
                if np.any(np.all(j_sort == np.sort(l), axis=1)):
                    incl.append(cl)
        vc_b += [len(np.unique(incl)) == num_dists]

    # eliminating false positives
    valid_combinations = valid_combinations[vc_b]

    if len(valid_combinations) > 0:
        # coords and data of valid motif
        coords_iv = res_coords[valid_combinations]
        data_iv = res_data[valid_combinations]
        if pymol_check:
            pymol_print(data_iv, coords_iv)
        if feedback:
            for ci, i in enumerate(data_iv):
                print(f"Motif {ci}")
                for k in i:
                    print("-".join(k))
                print()
        return data_iv, coords_iv

    else:
        if feedback:
            print("motif not present")
        return


if __name__ == "__main__":
    trianglex(
        "/PATH/TO/PDB/FILE",
        ["LEU", "ILE", "ILE", "LEU"],
        [6.11, 3.08, 5.85, 9, 9, 9],
        cushion=4,
        feedback=True,
    )
