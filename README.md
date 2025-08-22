
[User Manual](inst/docs/FUSION_v1.0.2_User_Manual.pdf)

# FUSION (Family-level Unique Small RNA Integration) Version 1.0.2

## Install instructions:

#### Directly from github

remotes::install_github("cozyrna/FUSION")

#### Using source code

1. Download FUSION_*.tar.gz file
2. Open Rstudio or R and type as below:

install.packages("~/FUSION_*.tar.gz", repos = NULL, type = "source")

library(FUSION)

#### There are three functions in this package : 

#### FUSION_ps - Diferential abundance analysis of sncRNA families in paired data samples using expression matrix

#### FUSION_ms - Diferential abundance analysis of sncRNA families in multiple samples data using expression matrix

#### FUSION_msmc - Diferential abundance analysis of sncRNA families in multiple samples data with multiple conditions using expression matrix
 
#### Ask help for the description and help menus of each of these functions :

?FUSION_ps

?FUSION_ms

?FUSION_msmc 

#### Note: If want to run examples straight in the Help documentation, it is necessary to first set working directory to the base folder of the installed package FUSION. Such as : setwd("/home/..../R/.../FUSION/"), or simply provide the exact path of the appropriate example matrix files.


#### Here are the sample runs with example data sets:

FUSION_ps(a = "./extdata/example_matrix_p1.txt", row_mean = 0.5)

FUSION_ms(a = "./extdata/example_matrix1.txt", S1 = 10, S2 = 16, row_mean = 1, top_species = 5000)

FUSION_msmc(a = "./extdata/example_matrix_cl.txt", cl = "./extdata/example_condition1.txt", row_mean = 1, top_species = 5000) 

#### The R script `plot_fusion.R` (available in the installed package under `scripts/`) provides a simple way to visualize the position of dysregulated RNA species along a parental RNA.
#### To use:
```r
# Load the script from the installed package
script_path <- system.file("scripts", "plot_fusion.R", package = "FUSION")

# Open and edit input paths if you want to use your own data
file.edit(script_path)  

# Run the script
source(script_path)

# Make sure to update in the script:
# - The path to your FASTA file
# - The path to your expression data file
# - The desired coordinates
# - The number of dysregulated RNA species to display (default: 100)

```
#### DOI: 10.5281/zenodo.16281885 (https://doi.org/10.5281/zenodo.16281885)
