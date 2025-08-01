% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FUSION.R
\name{FUSION_msmc}
\alias{FUSION_msmc}
\title{FUSION_msmc - Differential abundance analysis of sncRNA families in multiple samples data with multiple conditions using expression matrix}
\usage{
FUSION_msmc(
  a,
  cl,
  row_mean = getOption("row_mean", "0.1"),
  sncrna_family = c("tsrna", "rsrna", "ysrna", "mirna", "other"),
  top_species = getOption("max_num", "1000"),
  padj_method = c("bonferroni", "BH"),
  unique_anno = FALSE
)
}
\arguments{
\item{a}{a matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Rest of the columns provides RPM or other normalized expression values from different samples under study. Sequence (or ID) must be unique.}

\item{cl}{a file specifying multiple sample conditions in a comma separated format (such as for 3 different conditions with 6 samples each the file will contains the input: 1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3). Note: one can specify the condition in any order depending on the order of  the samples in the matrix, such as : “2,2,2,1,1,1,3,3,3” or “1,2,1,2,1,2,1,2,1,2” or “1,1,1,2,2,2,1,1,1,3,3,3”.}

\item{row_mean}{mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix}

\item{sncrna_family}{list of sncRNA families to be analyzed in the abundance analysis study. Use "tsrna", "rsrna", "ysrna", or "mirna" for tsRNAs, rsRNAs, ysRNAs, or mirna, respectively. Use "other" for 'pRNA,snRNA and snoRNA'. For all, use any letter or number, eg. "a","b", "c", 1, 2, 3. By default (i.e. if no option specified) it will search for tryRNAs (tsRNAs, rsRNAs, and ysRNAs)}

\item{top_species}{number (default 1000) of top species for each sncRNA family to be considered for analysis. It can help to reduce the run time of the analysis for the families (such as rsrna families) with large number of species. If time is not a concern, use high values such as 5000, 10000, etc.}

\item{padj_method}{adjustment methods for correcting p-values. Use either "bonferroni" or "BH", where BH stands for Benjamini & Hochberg. By default (i.e. if no option specified) it will run for "bonferroni".}

\item{unique_anno}{filter the input matrix to consider only uniquely mapped reads (sncRNA species), meaning only the sncRNA species that map to a single parent RNA. By default, it is FALSE (i.e., all relevant sncRNA species in the matrix are considered, irrespective of whether they map to a single or multiple parent RNAs).}
}
\value{
It will return a final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
}
\description{
FUSION_msmc - Differential abundance analysis of sncRNA families in multiple samples data with multiple conditions using expression matrix
}
\examples{
# Note:  After installation, one can find the example files in "../FUSION/extdata/"
# To run on your own matrix file, provide the full path as : 
# FUSION(a = "/path/to/your_matrix.txt", cl = "/path/to/your_condition.txt")
example_matrix_cl <- system.file("extdata", "example_matrix_cl.txt", package = "FUSION")
example_condition1 <- system.file("extdata", "example_condition1.txt", package = "FUSION")
example_condition2 <- system.file("extdata", "example_condition2.txt", package = "FUSION")
example_condition3 <- system.file("extdata", "example_condition3.txt", package = "FUSION")
# Run differential abundance analysis on example_matrix_cl.txt with 18 samples as per
# the conditions specified (i.e., 3 different conditions with 6 samples each) in the file 
# example_condition1.txt at row_mean threshold of 1 for 5000 top_species for 
# sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
FUSION_msmc(a = example_matrix_cl, cl = example_condition1, row_mean = 1, top_species = 5000) 
# Run differential abundance analysis on example_matrix_cl.txt with 18 samples as per
# the conditions specified (i.e., 4 different conditions with condition 1, 2, and 3 having five
# samples each, while last three samples are representing the condition 4) in the file 
# example_condition2.txt at row_mean threshold of 1 for 5000 top_species for 
# sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
FUSION_msmc(a = example_matrix_cl, cl = example_condition2, row_mean = 1, top_species = 5000)   
# Run differential abundance analysis on example_matrix_cl.txt with 18 samples as per
# the conditions specified (i.e., 2 different conditions with condition 1, and 2 having  10 and 
# 8 samples, respectively) in the file example_condition3.txt at row_mean threshold of 10 for 
# default (1000) top_species for tsRNA families. 
FUSION_msmc(a = example_matrix_cl, cl = example_condition3, row_mean = 10, sncrna_family = "tsrna") 
# Note:  If you want to save the terminal/console output to a file, use sink() command. 
# e.g.,:  options(max.print = 1e6); sink("~/output.txt"); 
# FUSION_ms(a = "./extdata/example_matrix1.txt", S1 = 10, S2 = 16); sink()
}
