<p align="center">
SCBD
</p>

## Description

This script computes the SCBD values for all features of a given normalised matrix, as described in <a href="https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.12141">Legendre and De Càceres 2013</a>, after tranformation of normalised counts using the Hellinger transformation.
The SCBD values give an estimate of the feature contributions to the overall signal, allowing to extract the most contributing features of a large dataset. The features can be species abundances, gene expression levels or any feature quantified by sequencing data. 


## Requirement
- Python version 3


## Usage
_ the script takes as input a <b>coma-separated</b> matrix of <b>normalised</b> counts:
```
python3 SCBD.py -m normalised_matrix.csv
```
It returns 2 output files:\n
\t_ out_SCBD_values.tsv\ttab-delimited file containing the SCBD values, sorted by decreasing order
\t_ out_SCBD_matrix.csv\tyour normalised matrix with features sorted by decreasing SCBD values\n
_ have a look at the distribution of SCBD values and determine the threshold below which the features should be discarded. It usually corresponds to the inflexion point of the distribution.\n
_ remove from the sorted matrix the features with SCBD values below the chosen threshold (the bottom part of the sorted matrix)
