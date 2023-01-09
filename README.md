![TEST](https://github.com/ugSUBMARINE/motif-search/actions/workflows/test.yml/badge.svg)

This repository contains a program to find a motif (eg active site constellation) of interest in a protein structure.
In order to find the motif, pseudoatom positions for each residues side chain are created. These represent the mean position of the catalytically important atoms of a side chain.
The input needs to be a pdb file of the protein.
![alt text](https://github.com/ugSUBMARINE/motif-search/blob/master/test_data/motif_sample.png?raw=true)


**Software Requirements:**
*  [Python3.10](https://www.python.org/downloads/)

*optional:*
*  [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/index.html)

In order to see all parameters run:

`python3 run_search.py -h`

Example usage:

`python3 run_search.py -f /PATH/TO/PDB/FILE -m SER-HIS-ASP -d 4.0-3.2-5.9` 
