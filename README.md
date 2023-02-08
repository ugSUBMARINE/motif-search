![TEST](https://github.com/ugSUBMARINE/motif-search/actions/workflows/test.yml/badge.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

This repository contains a program to find a motif (eg active site constellation) of interest in a protein structure.
In order to find the motif, pseudoatom positions for each residues side chain are created. These represent the mean position of the catalytically important atoms of a side chain.
The input needs to be a pdb file of the protein.
![alt text](https://github.com/ugSUBMARINE/motif-search/blob/master/test_data/motif_sample.png?raw=true)


**Software Requirements:**
*  [Python3.10](https://www.python.org/downloads/)

*optional:*
*  [Anaconda or Miniconda](https://docs.anaconda.com/anaconda/install/index.html)

In order to install the required packages run:
```
python3 -m venv /PATH/TO/YOUR/VENV
source venv/bin/activate
pip install -r requirements.txt
```

In order to see all parameters run:

`python3 run_search.py -h`

Example usage for a catalytic triad where the distance is 4.0A between SER-HIS, 3.2A between SER-ASP and 5.9A between HIS-ASP:

`python3 run_search.py -f /PATH/TO/PDB/FILE -m SER-HIS-ASP -d 4.0-3.2-5.9` 
