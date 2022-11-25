import urllib.request
from itertools import combinations

from bs4 import BeautifulSoup
import numpy as np

X_DICT = {
    "ALA": ["CB"],
    "CYS": ["SG"],
    "ASP": ["OD2", "OD1"],
    "GLU": ["OE1", "OE2"],
    "PHE": ["CG", "CD1", "CE1", "CD2", "CE2", "CZ"],
    "GLY": ["CA"],
    "HIS": ["CG", "ND1", "CE1", "NE2", "CD2"],
    "ILE": ["CB", "CG1", "CD1", "CG2"],
    "LYS": ["NZ"],
    "LEU": ["CB", "CG", "CD1", "CD2"],
    "MET": ["SD", "CE"],
    "ASN": ["OD1", "ND2"],
    "PRO": ["N", "CA", "CB", "CG", "CD"],
    "GLN": ["OE1", "NE2"],
    "ARG": ["NE", "CZ", "NH1", "NH2"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "VAL": ["CB", "CG1", "CG2"],
    "TRP": ["CE3", "CZ3", "CH2", "CZ2", "CE2", "NE1", "CD1", "CG", "CD2"],
    "TYR": ["OH"],
}


def access_url(url: str, parser: str):
    """access url of a website
    :parameter
        - url:
          url of interest
        - parser:
          which parser to use for the website eg 'html.parser'
    :return
        - content_parsed
          the websites content
    """
    # opens url of according uniprot entry
    acc_url = urllib.request.urlopen(url)
    # reads the data of the site
    content = acc_url.read()
    # parses the read data
    content_parsed = BeautifulSoup(content, parser)
    return content_parsed


def motif_dict_append(
    name: str, residues: list[str], distances: list[int | float]
) -> None:
    """appends defined motifs to the storage file
    :parameter
        - name:
          name of the motif
        - residues
          residues that make the motif
        - distances
          distances between all residues
    :return
        - None
    """
    with open("motif_dict.csv", "a") as motif_dict:
        motif_dict.write(
            name
            + ","
            + "_".join(residues)
            + ","
            + "_".join(np.asarray(distances, dtype=str))
            + "\n"
        )


def read_motif_dict() -> None:
    """reads in motif specifying file and converts it to dict
    :parameter
        - None:
    :return
        - m_dict:
          specifying the motifs - their residues and distances
    """
    m_dict = {}
    # reads all entires in the file and converts it to the needed name:[Res, dist]
    with open("motif_dict.csv", "r") as motif_dict:
        for ci, i in enumerate(motif_dict):
            if ci > 0:
                i = i.strip()
                i_split = i.split(",")
                m_dict[i_split[0]] = [
                    i_split[1].split("_"),
                    np.asarray(i_split[2].split("_"), dtype=float),
                ]
    return m_dict


def pseudoatom_positions(
    atom_pos_dict: dict, target_pdb_file: str | None = None, pdb_code: str | None = None
) -> tuple[
    np.ndarray[tuple[int, int], np.dtype[str]],
    np.ndarray[tuple[int, int], np.dtype[float]],
]:
    """calculates pseudoatom positions for all residues in target_pdb_file
        which is the mean distance between all atoms according to atom_pos_dict
     :parameter
        - atom_pos_dict: every amino acid with its catalytically important atoms eg
          {"ASP": ["OD2", "OD1"],...}
        - target_pdb_file:
          pdb file with data of protein of interest
        - pdb_code:
          pdb code for data of protein of interest from the pdb website
          without saving the file and target_pdb_file will be ignored
    :return
        pseudo_data: 2D list like [[Res 3letter, ChainID, ResidueID],...]
        pseudo_coords: 2D list like [[pseudo_x1, pseudo_y1, pseudo_z1],...]"""
    # list of all data of the entries like
    # [[Atom type, Residue 3letter, ChainID, ResidueID],...]
    res_data = []
    # list of all coordinates of the entries like [[x1, y1, z1],...]
    res_coords = []
    if pdb_code is None:
        # read all lines
        file = open(target_pdb_file, "r")
    else:
        # parse information from website
        web_acc = access_url(
            "https://files.rcsb.org/view/" + pdb_code.lower() + ".pdb",
            "html.parser",
        )
        # format the data so it can be used in the for loop
        file = str(web_acc).split("\n")
    for line in file:
        if "ATOM  " in line[:6]:
            line = line.strip()
            res_data += [
                [
                    line[12:16].replace(" ", ""),
                    line[17:20].replace(" ", ""),
                    line[21].replace(" ", ""),
                    line[22:26].replace(" ", ""),
                ]
            ]
            res_coords += [[line[30:38], line[38:46], line[46:54]]]
    if pdb_code is None:
        file.close()

    res_data = np.asarray(res_data)
    res_coords = np.asarray(res_coords, dtype=float)

    aa = [
        "GLY",
        "ALA",
        "VAL",
        "LEU",
        "ILE",
        "PRO",
        "PHE",
        "TYR",
        "TRP",
        "SER",
        "THR",
        "MET",
        "CYS",
        "ASP",
        "ASN",
        "GLU",
        "GLN",
        "ARG",
        "LYS",
        "HIS",
    ]
    # list of 2D arrays where each 2D array is like
    # [[Res 3letter, ChainId, ResidueID],...]
    # each 2D array is for an aa
    pseudo_data = []
    # list of 2D arrays where each 2D array is like
    # [[pseudo_x1, pseudo_y1, pseudo_z1],...]
    # each 2D array is for an aa
    pseudo_coords = []
    for i in aa:
        # where aa i is located i res_data and res_coords
        ind_in_ori_arr = np.where(res_data[:, 1] == i)[0]
        # if aa i exists
        if len(ind_in_ori_arr) > 0:
            # which atom is a catalytically important atom
            cat_imp_at_ind = np.isin(res_data[ind_in_ori_arr][:, 0], atom_pos_dict[i])

            # coordinates[where aa in res_data/coords][catalytically important atom]
            cat_imp_at_coords = res_coords[ind_in_ori_arr][cat_imp_at_ind]
            cat_im_at_data = res_data[ind_in_ori_arr][cat_imp_at_ind]

            # to get for each residue one entry with [Res 3letter, ChainID, ResID]
            pseudo_data += np.asarray(
                np.split(cat_im_at_data, len(cat_imp_at_coords) / len(atom_pos_dict[i]))
            )[:, :, 1:4][:, 0].tolist()

            # pseudo coordinates for all residues of aa i
            pseudo_coords += np.mean(
                np.asarray(
                    np.split(
                        cat_imp_at_coords,
                        len(cat_imp_at_coords) / len(atom_pos_dict[i]),
                    )
                ),
                axis=1,
            ).tolist()

    return pseudo_data, pseudo_coords


def motif_extraction(
    res_oi_ind: list[int],
    res_oi_chain: list[str],
    motive_template_path: None | str = None,
    pdb_code: None = None,
    precision: int = 2,
):
    """extracts the three letter code of the res_oi_nid and the distances of each
    residues pseudoatom position against each other
    :parameter
        - res_oi_ind:
          indices of the residues that form the motive eg [24, 50, 70]
        - res_oi_chain:
          chain ids of the residues that form the motive eg ["A", "A", "A"]
        - motive_template_path:
          path to the pdb file of the protein which serves as a template to get the
          distances between the residues that form the motive of interest
        - pdb_code:
          pdb code for data of protein of interest from the pdb website
          without saving the file and target_pdb_file will be ignored
        - precision
          number of decimal places of the measured distances
    :return
        - three_letter:
          three letter code of the residues of interest eg ['ASP', 'GLU', 'VAL']
        - comb_dist_between_res:
          distances between the residues of interest
          eg [38.44, 26.00, 33.31]
          distances are like
          from itertools import combinations
          list(combinations(range(len(res_oi_ind)), 2)) which results in
          [(0, 1), (0, 2), (1, 2)] what means the first distance is between ASP-GLU,
          second between ASP-VAL and the third between GLU-VAL
    """

    if len(res_oi_ind) != len(res_oi_chain):
        print(
            "Error residue IDs and chain IDs don't match\ngiven residue IDs: "
            + str(len(res_oi_ind))
            + " - but given chain IDs: "
            + str(len(res_oi_chain))
        )
        return
    pdb_data = np.asarray(
        pseudoatom_positions(X_DICT, motive_template_path, pdb_code=pdb_code)
    )
    # all possible pairs between the res_oi
    res_combs = list(combinations(range(len(res_oi_ind)), 2))

    # the coordinates and the three letter code of the residues of interest (res_oi)
    pseudo_coords = []
    three_letter = []
    for i, j in zip(res_oi_ind, res_oi_chain):
        # where the res_oi with the right ind and the right chain combined in located
        needed_bool_id = np.isin(np.asarray(pdb_data[0][:, 2], dtype=int), i)
        needed_bool_chain = np.isin(np.asarray(pdb_data[0][:, 1], dtype=str), j)
        needed_bool = np.all((needed_bool_chain, needed_bool_id), axis=0)
        if len(pdb_data[0][needed_bool][:, 0].tolist()) == 0:
            raise ValueError(
                f"Residue 'chain {j} index {i}' is not present in this structure"
            )
            return
        three_letter += pdb_data[0][needed_bool][:, 0].tolist()
        pseudo_coords += pdb_data[1][needed_bool].tolist()

    pseudo_coords = np.asarray(pseudo_coords, dtype=float)

    # get only the x,y,z coordinates from the input arrays and reshape them so
    # they can be subtracted from each other
    ind_coord_arr1_rs = pseudo_coords.reshape(pseudo_coords.shape[0], 1, 3)
    ind_coord_arr2_rs = pseudo_coords.reshape(1, pseudo_coords.shape[0], 3)
    # calculating the distance between each point
    distance = np.round(
        np.sqrt(np.sum((ind_coord_arr1_rs - ind_coord_arr2_rs) ** 2, axis=2)), precision
    )

    # distance according to the formed pairs in res_combs
    comb_dist_between_res = []
    for i in res_combs:
        comb_dist_between_res += [distance[i[0]][i[1]]]

    return [three_letter, comb_dist_between_res]


def pymol_print(data, coords):
    intermediate_sticks = []
    # show sticks of the valid_combinations residues
    for i in np.unique(data.reshape((-1, 3))[:, [1, 2]], axis=0):
        intermediate_sticks += ["".join(["(chain ", i[0], " and resi ", i[1], ")"])]
    print("show sticks, ", " or ".join(intermediate_sticks))
    # show their pseudoatom position
    for i in np.unique(coords.reshape((-1, 3)), axis=0):
        print(
            "pseudoatom tmpPoint2, resi=40, chain=ZZ, b=40, color=tv_blue, pos=",
            i.tolist(),
        )


if __name__ == "__main__":
    pass
    # print(motif_extraction([4, 6, 15], ["A", "A", "A"], "/home/gwirn/gb1.pdb"))
    pseudoatom_positions(X_DICT, pdb_code="4JYM")
