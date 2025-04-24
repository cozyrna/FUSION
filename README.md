# FUSION (Family-level Unique Small RNA Integration) Version 1.0.0

## Install instructions:

#### Directly from github

remotes::install_github("cozyrna/FUSION")

#### Using source code

1. Download FUSION_*.tar.gz file
2. Open Rstudio or R and type as below:

install.packages("~/FUSION_*.tar.gz", repos = NULL, type = "source")

library(FUSION)

#### There are three functions in this package : 

#### FUSION_ps - Enrichment Analysis of sncRNA families in paired data samples using expression matrix

#### FUSION_ms - Enrichment Analysis of sncRNA families in multiple samples data using expression matrix

#### FUSION_msmc - Enrichment Analysis of sncRNA families in multiple samples data with multiple conditions using expression matrix
 
#### Ask help for the description and help menus of each of these functions :

?FUSION_ps

?FUSION_ms

?FUSION_msmc 

#### Note: If want to run examples straight in the Help documentation, it is necessary to first set working directory to the base folder of the installed package FUSION. Such as : setwd("/home/..../R/.../FUSION/"), or simply provide the exact path of the appropriate example matrix files.

setwd("~/R/x86_64-pc-linux-gnu-library/4.3/FUSION")

### Here are the sample runs with example data sets:

FUSION_ps(a = "./extdata/example_matrix_p1.txt", row_mean = 0.5)

FUSION_ms(a = "./extdata/example_matrix1.txt", S1 = 10, S2 = 16, row_mean = 1, top_species = 5000)

FUSION_msmc(a = "./extdata/example_matrix_cl.txt", cl = "./extdata/example_condition1.txt", row_mean = 1, top_species = 5000) 

