# R scripts to perform single-cell eQTL mapping in Boocock et al. 2024

## Analysis of segregants

The primary script for analyzing the segregant data, including genotyping, eQTL mapping, identifying hotspots, mapping cell-cycle occupancy, and identifying cell-cycle X eQTL interactions can be found here: [processSegregants.R](processSegregants.R)  

The R variable `code.dir` should be set to a folder containing the R scripts downloaded from [yeast/code](https://github.com/joshsbloom/single_cell_eQTL/tree/master/yeast/code)  

Data in processed/ reference/ and cell cycle annotations can be found on [Dryad](https://doi.org/10.5061/dryad.xgxd254qb)  

## Analysis of hybrid diploids

The primary script for calculating allele-specific expression mean and dispersion effects can be found here: [processDiploids.R](processDiploids.R)

Allele-specific expression dispersion power calculation code can be found here: [simASEdisp.R](simASEdisp.R)
