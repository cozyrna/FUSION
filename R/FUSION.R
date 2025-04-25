#!/usr/bin/Rscript

##' @name FUSION
##' @rdname FUSION
##' @title FUSION  - Family-level Unique Small RNA Integration across samples using expression matrix
##' @details
##' For more information, see the user manual:
##' \code{browseURL(system.file("docs/FUSION_1.0_User_Manual.pdf", package = "FUSION"))}
##' 
##' After installing the package FUSION, call the library as: library(FUSION).
##' \cr
##' \cr There are three main functions in this package : FUSION_ps, FSION_ms, and FUSION_msmc
##' \cr
##' \cr Ask help for the description and help menus of each of these functions :
##' \cr ?FUSION_ps
##' \cr ?FUSION_ms
##' \cr ?FUSION_msmc
##' \cr
##' \cr If want to run examples, it is necessary to set working directory to the base folder of the installed package FUSION. Such as : setwd("/home/..../R/.../FUSION/").
##' \cr
##' @param a a matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Sequence (or ID) must be unique.  Rest of the columns provides RPM or expression values from different samples under study
##' @param row_mean  mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
##' @param sncrna_family list of sncRNA families to be analyzed for expression analysis study. Use "tsrna", "rsrna", "ysrna", or "mirna" for tsRNAs, rsRNAs, ysRNAs, or mirna, respectively. Use "other" for 'pRNA, snRNA and snoRNA'. For all, use any letter or number, eg. "a","b", "c", 1, 2, 3. By default (i.e. if no option specified) it will search for tryRNAs (tsRNAs, rsRNAs, and ysRNAs)
##'
##' @return with FUSION_ps, it will return, for each pair in the input matrix, an output data-frame with w_positive, w_negative, P_value, and adjusted P_value for each sncRNA family chosen for analysis.
##' @return with FUSION_ms, it will return a final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
##' @return with FUSION_msmc, it will return a final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
#'
#' @importFrom utils read.delim
#' @importFrom usethis use_package
#' @importFrom stats lm p.adjust wilcox.test

NULL

#' FUSION_ps - Differential expression analysis of sncRNA families in paired data samples using expression matrix
#'
#' @param a a matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Sequence (or ID) must be unique.  Rest of the columns provides RPM or expression values from different samples under study, such as first half of those belongs to the samples from Condition1 (such as control or from healthy tissue) and second half belongs to the samples from Condition2 (such as treated or infected tissue)
#' @param order use either G or P for specifying the order of paired samples in the input matrix, where G is for Samples in Group (i.e. all samples from Condition1 followed by all samples from Condition2) (follow example_matrix_p1)) and P is for samples in Pairs (follow example_matrix_p4). By default it considers sample pairs in Group (G)")
#' @param row_mean  mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
#' @param sncrna_family list of sncRNA families to be analysis for expression analysis study. Use "tsrna", "rsrna", "ysrna", or "mirna" for tsRNAs, rsRNAs, ysRNAs, or mirna, respectively. Use "other" for 'pRNA,snRNA and snoRNA'. For all, use any letter or number, eg. "a","b", "c", 1, 2, 3. By default (i.e. if no option specified) it will search for tryRNAs (tsRNAs, rsRNAs, and ysRNAs)
#' @param padj_method adjustment methods for correcting p-values. Use either "bonferroni" or "BH", where BH stands for Benjamini & Hochberg. By default (i.e. if no option specified) it will run for "bonferroni". 
#' @return For each pair in the input matrix, it will return an output data-frame with w_positive, w_negative, P_value, and adjusted P_value for each sncRNA family chosen for analysis.
#'
#' @examples
#' # Note:  After installation, one can find the example files in "../FUSION/extdata/"
#' # To run on your own matrix file, provide the full path as : FUSION(a = "/path/to/your_matrix.txt")
#' example_matrix_p1 <- system.file("extdata", "example_matrix_p1.txt", package = "FUSION")
#' example_matrix_p2 <- system.file("extdata", "example_matrix_p2.txt", package = "FUSION")
#' example_matrix_p3 <- system.file("extdata", "example_matrix_p3.txt", package = "FUSION")
#' example_matrix_p4 <- system.file("extdata", "example_matrix_p4.txt", package = "FUSION")
#' # Run differential expression analysis on example_matrix_p1.txt (5 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
#' FUSION_ps(a = example_matrix_p1) 
#' # Run differential expression analysis on example_matrix_p1.txt (5 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
#' FUSION_ps(a = example_matrix_p1) 
#' # Run differential expression analysis on example_matrix_p1.txt (5 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default) using
#' # BH (Benjamini & Hochberg) method for correcting or adjusting p-values;
#' FUSION_ps(a = example_matrix_p1, padj_method = "BH") 
#' # Run differential expression analysis on example_matrix_p1.txt (5 pairs of samples) at row_mean
#' # threshold of 0.5 for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
#' FUSION_ps(a = example_matrix_p1, row_mean = 0.5) 
#' # Run differential expression analysis on example_matrix_p2.txt (10 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for all sncRNA families;
#' FUSION_ps(a = example_matrix_p2, sncrna_family = "a")
#' # Run differential expression analysis on example_matrix_p2.txt (10 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for all sncRNA families;
#' FUSION_ps(a = example_matrix_p2, sncrna_family = 0)
#' # Run differential expression analysis on example_matrix_p2.txt (10 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for miRNA families;
#' FUSION_ps(a = example_matrix_p2, sncrna_family = "mirna")
#' # Run differential expression analysis on example_matrix_p2.txt (10 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for rsRNA families;
#' FUSION_ps(a = example_matrix_p2, sncrna_family = "rsrna") 
#' # Run differential expression analysis on example_matrix_p2.txt (4 pairs of samples) at row_mean
#' # threshold of 0.5 for tsRNA families;
#' FUSION_ps(a = example_matrix_p3, row_mean = 0.5, sncrna_family = "tsrna")  
#' # Run differential expression analysis on example_matrix_p2.txt (4 pairs of samples) at row_mean
#' # threshold of 0.1 for for ysRNA families;
#' FUSION_ps(a = example_matrix_p3, row_mean = 0.1, sncrna_family = "ysrna")  
#' # Run differential expression analysis on example_matrix_p2.txt (4 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for other (pRNA,snRNA and snoRNA) sncRNA families;
#' FUSION_ps(a = example_matrix_p3, sncrna_family = "other")  
#' # Run differential expression analysis on example_matrix_p4.txt (5 pairs of samples) at default
#' # row_mean threshold (i.e., 0.1) for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default) and 
#' # samples are arranged in as pairs of columns;
#' FUSION_ps(a = example_matrix_p4,order = "P") 
#' # Note:  If you want to save the terminal/console output to a file, use sink() command. e.g.:
#' # options(max.print = 1e6); sink("~/output.txt"); FUSION_ps(a = "example_matrix_p1.txt"); sink()
#' @export

  FUSION_ps <- function(a, row_mean = getOption("row_mean", "0.1"), sncrna_family = c("tsrna", "rsrna", "ysrna", "mirna", "other"), padj_method = c("bonferroni", "BH"), order = c("G", "P")) {
  #### FUSION_ps function body starts ####
  args = commandArgs(trailingOnly=TRUE)
  
  ########## Part1: defining the function to call wilcox.test for a matrix of a pair of samples ##########
  
  options(warn=1)
  sncrna.family.paired_wilcox <- function(e1, e2)
  {
    diff = e1 - e2
    diff = diff[diff != 0]
    r = rank(abs(diff))
    signed_r = r * sign(diff)
    w_pos = sum(signed_r[signed_r>0])
    w_neg = abs(sum(signed_r[signed_r<0]))
    p = wilcox.test(e1, e2, paired=T, exact = FALSE)$p.value
    return(c(w_pos, w_neg, p))
  }
  
  ########## Part2: Processing of input matrix and different parameters to be used in the function ##########
  
  ########## 2.1: Preprocessing matrix if samples are arranged as pairwise columns ##########
  
  # A Function to rearrange columns dynamically by grouping
  rearrange_columns <- function(df) {
    # Identify the total number of columns
    num_cols <- ncol(df)
    
    # Assuming columns are in pairs and alternate between two categories
    # For example, first two columns are Sequence and annotation, while columns 3 and 4 makes a pair, columns 5 and 6 makes another pair, so on.
    
    # Split columns into pairs
    
    # Columns for the first category
    group_1 <- seq(3, num_cols, by = 2)  # Even columns (3, 5, 7, ...)
    
    # Columns for the second category
    group_2 <- seq(4, num_cols, by = 2)  # Odd columns (4, 6, 8, ...)
    
    # Reorder columns: first the first 2 columns, then the first category, and finally the second category
    df_rearranged <- df[, c(1:2, group_1, group_2)]
    
    return(df_rearranged)
  }
  
  #### a: a user-provided matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Rest of the columns provides RPM or expression values from different samples under study. Sequence (or ID) must be unique.####
  a <- read.delim(a, header = TRUE)
  a <- as.data.frame(a)
  
  if (missing(order)) {
    #order = "G"
    a <- a
  } else if(order == "G"){ 
    #order = "G"
    a <- a
  } else if(order == "P"){ 
    #order = "P"
    # Apply the rearrange_columns function
    a <- rearrange_columns(a)
  } else{
    print("Error: use either G or P for specifying the order of paired samples in the input matrix, where G is for Samples in Group (i.e. all samples from group 1 followed by all samples from group 2) and P is for samples are in Pairs. By default it considers pairs in groups (see example matrix)")
    return(NULL)
  }
  
  #### e and anno are sub-matrix generated from the user-provided input matrix (a)
  #### e : An expression matrix with first column as Sequence (or ID) and rest of the columns with expression values
  #### anno :  Sequence (or ID) and their Annotation information
  
  anno = a[,1:2]  #### Separated the columns with unique Sequence (or ID) and their Annotation information
  e = a[,3:ncol(a)] #### Separated the columns with expression values
  rownames(e) = a[,1]
  
  if (ncol(e) %% 2 != 0){
    print("Error: column numbers are not in pair in the matrix. Please check")
    return(NULL)
  }
  
  ##### Extracting unique annotation list from input file
  uni_annotation = unique(unlist(strsplit(format(anno[,2], justify="none"), ";")))
  
  #### Extracting genomic tRNAs in the annotation list
  tmp1 = uni_annotation[grepl("mature-tRNA", uni_annotation)]
  gtsrna_family = sort(unique(substr(tmp1, 1, 19)))
  gtsrna_family[gtsrna_family=="mature-tRNA-iMet-CA"] = "mature-tRNA-iMet-CAT"
  
  #### Extracting mitochondrial tRNAs in the annotation list
  tmp2 = uni_annotation[grepl("mature-mt_tRNA", uni_annotation)]
  mtsrna_family = sort(unique(substr(tmp2, 1, 22)))
  
  #### Extracting rRNAs in the annotation list
  #rsrna_family = c("4.5S-rRNA", "5S-rRNA", "5.8S-rRNA", "12S-rRNA", "16S-rRNA", "18S-rRNA", "28S-rRNA", "45S-rRNA")
  tmp3 = uni_annotation[grepl("S-rRNA", uni_annotation)]
  rsrna_family = sort(unique(substr(tmp3, 1, 9)))
  
  #### Extracting YRNAs in the annotation list
  #ysrna_family = c("RNY1-YRNA", "RNY3-YRNA", "RNY4-YRNA", "RNY5-YRNA")
  tmp4 = uni_annotation[grepl("-YRNA", uni_annotation)]
  ysrna_family = sort(unique(substr(tmp4, 1, 9)))
  
  #### Extracting miRNAs in the annotation list
  mirna_patterns <- c("-lin", "-let", "-mir", "-miR", "-MIR")
  tmp5 = uni_annotation[grepl(paste(mirna_patterns, collapse='|'), uni_annotation)]
  mirna_family = sort(unique(substr(tmp5, 1, 16)))
  
  #### Extracting other sncRNAs (pRNA,snRNA and snoRNA) in the annotation list
  other_family <- c("piRNA", "snRNA", "snoRNA")
  
  if (missing(sncrna_family)) {
    sncrna_family = c(gtsrna_family, mtsrna_family, rsrna_family, ysrna_family)
  } else if(sncrna_family == "tsrna"){
    sncrna_family = c(gtsrna_family, mtsrna_family)
  } else if(sncrna_family == "rsrna"){
    sncrna_family = rsrna_family
  } else if(sncrna_family == "ysrna"){
    sncrna_family = ysrna_family
  } else if (sncrna_family == "mirna"){
    sncrna_family = mirna_family
  } else if (sncrna_family == "other"){
    sncrna_family = other_family
  } else{
    sncrna_family = c(gtsrna_family, mtsrna_family, rsrna_family, ysrna_family, mirna_family, other_family)
  }
  
  #### row_mean : mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
  row_mean = as.integer(row_mean)
  
  e = e[rowMeans(e) > row_mean,]   #### only retain the sncRNA species with mean RPM > row_mean
  
  #### p_correction_method : adjustment method ( default : the Bonferroni correction) for correcting or adjusting the p-values
  #p_correction_method = "bonferroni"
  #padj_method = as.character(padj_method)
  if (missing(padj_method)) {
    padj_method = "bonferroni"
  } else if(padj_method == "bonferroni"){ 
    padj_method = "bonferroni"
  } else if(padj_method == "BH"){ 
    padj_method = "BH"
  } else{
    print("Error: use either bonferroni or BH for correcting or adjusting the p-values.")
    return(NULL)
  }
  
  ########## Part3: calling function to perform differential expression analysis on each pair using processed inputs and return output ##########
  
  ##### enlisting the sncRNA families (from Annotation column) in the matrix
  
  sncrna_family_new = c()
  
  if (length(sncrna_family) < 1) { 
    print ("There are no species annotated for the chosen sncRNA family in the input matrix at the given threshold")
    return(NULL)
  }
  else {
  
  
  for (k in 1:length(sncrna_family))
  {
    sncrna_family_c = grepl(sncrna_family[k], anno$Annotation) 
    if (grepl("rRNA", sncrna_family[k])) sncrna_family_c = (anno$Annotation==sncrna_family[k])
    
    else if(grepl("mir", sncrna_family[k], ignore.case = TRUE)) sncrna_family_c = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
    else if(grepl("let", sncrna_family[k])) sncrna_family_c = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
    
    else sncrna_family_c = grepl(sncrna_family[k], anno$Annotation)
    e_tmp = e[sncrna_family_c,]
    count <- nrow(e_tmp)
    if(count==0){
      sncrna_family_new <- c(sncrna_family_new, sncrna_family[(k)])
    }
  }
  sncrna_family <- sncrna_family[!sncrna_family %in% sncrna_family_new]
  }
  ##### running the defined function for each pair in the matrix
  
  n = ncol(e)/2
  
  i <- 1
  output = c()
  while (i <= n){
    e_ij = data.frame(e[,i], e[, i+n])
    rownames(e_ij) = rownames(e)
    
    if (ncol(e_ij) != 2)
    {
      print("Error: column number is not muliple of two")
      return(NULL)
    }
    
    z = match(rownames(e_ij), anno$Sequence)
    if (length(z[is.na(z)]) > 0)
    {
      print("Error: annotation cannot match the input data")
      return(NULL)
    }
    anno = anno[z,]
    
    w_pos = rep(NA, length(sncrna_family))
    w_neg = rep(NA, length(sncrna_family))
    p = rep(NA, length(sncrna_family))
    
    for (k in 1:length(sncrna_family)){    
      if (grepl("rRNA", sncrna_family[k])) is.sncrna_family = (anno$Annotation==sncrna_family[k])
      
      else if(grepl("mir", sncrna_family[k], ignore.case = TRUE)) is.sncrna_family = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
      else if(grepl("let", sncrna_family[k])) is.sncrna_family = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
      
      else is.sncrna_family = grepl(sncrna_family[k], anno$Annotation)
      e_tmp = e_ij[is.sncrna_family,]
      out_tmp = sncrna.family.paired_wilcox(e_tmp[,1], e_tmp[,2])
      w_pos[k] = out_tmp[1]
      w_neg[k] = out_tmp[2]
      p[k] = out_tmp[3]
    }
    adjusted_p = p.adjust(p, method=padj_method)    #### Correcting or adjusting the p-values with padj_method (i.e. the Bonferroni correction (default) or by BH (Benjamini & Hochberg))
    output_tmp = data.frame(sncrna_family, w_pos, w_neg, p, adjusted_p)   #### writing final output in a data-frame with w_positive, w_negative, P_value, and adjusted P_value for each sncRNA family chosen for analysis
    output = rbind(output, data.frame(Pair = paste0("Pair_", i), output_tmp))
    i <- i+1
  }
  return(output)
}

########## FUSION_ps function ends ##########


#' FUSION_ms - Differential expression analysis of sncRNA families in multiple samples data using expression matrix
#'
#' @param a a matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Rest of the columns provides RPM or expression values from different samples under study. Sequence (or ID) must be unique.
#' @param S1 Number of samples from Condition1 (such as control or from healthy tissue)
#' @param S2 Number of samples from Condition2 (such as treated or infected tissue)
#' @param row_mean  mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
#' @param sncrna_family list of sncRNA families to be analysis for expression analysis study. Use "tsrna", "rsrna", "ysrna", or "mirna" for tsRNAs, rsRNAs, ysRNAs, or mirna, respectively. Use "other" for 'pRNA,snRNA and snoRNA'. For all, use any letter or number, eg. "a","b", "c", 1, 2, 3. By default (i.e. if no option specified) it will search for tryRNAs (tsRNAs, rsRNAs, and ysRNAs)
#' @param top_species number (default 1000) of top species for each sncRNA family to be considered for analysis. It can help to reduce the run time of the analysis for the families (such as rsrna families) with large number of species. If time is not a concern, use high values such as 5000, 10000, etc.
#' @param padj_method adjustment methods for correcting p-values. Use either "bonferroni" or "BH", where BH stands for Benjamini & Hochberg. By default (i.e. if no option specified) it will run for "bonferroni". 
#' @return It will return a final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
#'
#' @examples
#' # Note:  After installation, one can find the example files in "../FUSION/extdata/"
#' # To run on your own matrix file, provide the full path as : FUSION(a = "/path/to/your_matrix.txt")
#' example_matrix1 <- system.file("extdata", "example_matrix1.txt", package = "FUSION")
#' example_matrix2 <- system.file("extdata", "example_matrix2.txt", package = "FUSION")
#' example_matrix3 <- system.file("extdata", "example_matrix3.txt", package = "FUSION")
#' FUSION_ms(a = example_matrix1, S1 = 10, S2 = 16, row_mean = 1, top_species = 5000)   
#' # For running differential expression analysis on example_matrix1.txt with 10 healthy samples (S1) 
#' # and 16 patients (S2) at row_mean threshold of 1 for 5000 top_species for sncRNA families (tsRNAs,
#' # rsRNAs and ysRNAs) (default);
#' FUSION_ms(a = example_matrix1, S1 = 10, S2 = 16)    
#' # For running differential expression analysis on example_matrix1.txt with 10 samples from 
#' # Condition1 (S1) and 16 samples from Condition2 (S2) at default row_mean (i.e., 0.1) and 
#' # top_species threshold for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
#' FUSION_ms(a = example_matrix1, S1 = 10, S2 = 16, padj_method = "BH") 
#' # For running differential expression analysis on example_matrix1.txt with 10 samples from
#' # Condition1 (S1) and 16 samples from Condition2 (S2) at default row_mean (i.e., 0.1) and
#' # top_species threshold for sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default)
#' # using BH (Benjamini & Hochberg) method for correcting or adjusting p-values;
#' FUSION_ms(a = example_matrix1, S1 = 10, S2 = 16, sncrna_family = "a")    
#' # For running differential expression analysis on example_matrix1.txt with 10 samples from
#' # Condition1 (S1) and 16 samples from Condition2 (S2) at default row_mean threshold (i.e., 0.1)
#' # and top_species for all sncRNA families;
#' FUSION_ms(a = example_matrix2, S1 = 5, S2 = 5, sncrna_family = 0)    
#' # For running differential expression analysis on example_matrix2.txt with 5 healthy 
#' # samples (S1) and 5 patients (S2) at default row_mean threshold  (i.e., 0.1) and 
#' # top_species for all sncRNA families;
#' FUSION_ms(a = example_matrix2, S1 = 5, S2 = 5, sncrna_family = "mirna", top_species = 1000) 
#' # For running differential expression analysis on example_matrix2.txt with 5 samples from 
#' # Condition1 (S1) and 5 samples from Condition2 (S2) at default row_mean threshold (i.e., 0.1) 
#' # for 1000 top_species for miRNA families;
#' FUSION_ms(a = example_matrix2, S1 = 5, S2 = 5, row_mean = 10, top_species = 2000, 
#' sncrna_family = "rsrna")   
#' # For running differential expression analysis on example_matrix2.txt with 5 samples from 
#' # Condition1 (S1) and 5 samples from Condition2 (S2) at row_mean threshold of 10 for 
#' # 2000 top_species for rsRNA families;
#' FUSION_ms(a = example_matrix3, S1 = 10, S2 = 8, row_mean = 10, sncrna_family = "tsrna")   
#' # For running differential expression analysis on example_matrix3.txt with 10 control 
#' # samples (S1) and 8 treated samples (S2) at row_mean threshold of 10 for default (1000) 
#' # top_species for tsRNA families;
#' FUSION_ms(a = example_matrix3, S1 = 10, S2 = 8, row_mean = 0.1, sncrna_family = "ysrna")    
#' # For running differential expression analysis on example_matrix3.txt with 10 samples from 
#' # Condition1 (S1) and 8 samples from Condition2 (S2) at row_mean threshold of 0.1 for 
#' # default (1000) top_species for ysRNA families;
#' FUSION_ms(a = example_matrix3, S1 = 10, S2 = 8, top_species = 100, sncrna_family = "other")   
#' # For running differential expression analysis on example_matrix3.txt with 10 samples from 
#' # Condition1 (S1) and 8 samples from Condition2 (S2) at default row_mean threshold  (i.e., 0.1)
#' # for 100 top_species for other (pRNA,snRNA and snoRNA) sncRNA families;
#' # Note:  If you want to save the terminal/console output to a file, use sink() command. e.g.:  
#' # options(max.print = 1e6); sink("~/output.txt"); 
#' # FUSION_ms(a = "./extdata/example_matrix1.txt", S1 = 10, S2 = 16); sink()
#' @export

  FUSION_ms <- function(a, S1, S2, row_mean = getOption("row_mean", "0.1"), sncrna_family = c("tsrna", "rsrna", "ysrna", "mirna", "other"), top_species = getOption("max_num", "1000"), padj_method = c("bonferroni", "BH")){
  #### FUSION_ms function body starts ####
  args = commandArgs(trailingOnly=TRUE)
  
  
  ########## Part1: Processing of input matrix and different parameters to be used in the function ##########
  #### a: a user-provided matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Rest of the columns provides RPM or expression values from different samples under study. Sequence (or ID) must be unique.####
  a <- read.delim(a, header = TRUE)
  a <- as.data.frame(a)
  
  #### e and anno are sub-matrix generated from the user-provided input matrix (a)
  #### e : An expression matrix with first column as Sequence (or ID) and rest of the columns with expression values
  #### anno :  Sequence (or ID) and their Annotation information
  
  anno = a[,1:2]  #### Separated the columns with unique Sequence (or ID) and their Annotation information
  e = a[,3:ncol(a)] #### Separated the columns with expression values
  rownames(e) = a[,1]
  
  #### S1 and S2: Number of samples from Condition1 Condition2
  #### S1 : number of samples from Condition1
  S1 = as.numeric(S1)
  #### S2 : number of samples from Condition2
  S2 = as.numeric(S2)
  
  #### cl: a vector specifying samples from two conditions (1 and 2)
  cl <-c(rep(1, S1), rep(2, S2))
  
  ##### Extracting unique annotation list from input file
  uni_annotation = unique(unlist(strsplit(format(anno[,2], justify="none"), ";")))
  
  #### Extracting genomic tRNAs in the annotation list
  tmp1 = uni_annotation[grepl("mature-tRNA", uni_annotation)]
  gtsrna_family = sort(unique(substr(tmp1, 1, 19)))
  gtsrna_family[gtsrna_family=="mature-tRNA-iMet-CA"] = "mature-tRNA-iMet-CAT"
  
  #### Extracting mitochondrial tRNAs in the annotation list
  tmp2 = uni_annotation[grepl("mature-mt_tRNA", uni_annotation)]
  mtsrna_family = sort(unique(substr(tmp2, 1, 22)))
  
  #### Extracting rRNAs in the annotation list
  #rsrna_family = c("4.5S-rRNA", "5S-rRNA", "5.8S-rRNA", "12S-rRNA", "16S-rRNA", "18S-rRNA", "28S-rRNA", "45S-rRNA")
  tmp3 = uni_annotation[grepl("S-rRNA", uni_annotation)]
  rsrna_family = sort(unique(substr(tmp3, 1, 9)))
  
  #### Extracting YRNAs in the annotation list
  #ysrna_family = c("RNY1-YRNA", "RNY3-YRNA", "RNY4-YRNA", "RNY5-YRNA")
  tmp4 = uni_annotation[grepl("-YRNA", uni_annotation)]
  ysrna_family = sort(unique(substr(tmp4, 1, 9)))
  
  #### Extracting miRNAs in the annotation list
  mirna_patterns <- c("-lin", "-let", "-mir", "-miR", "-MIR")
  tmp5 = uni_annotation[grepl(paste(mirna_patterns, collapse='|'), uni_annotation)]
  mirna_family = sort(unique(substr(tmp5, 1, 16)))
  
  #### Extracting other sncRNAs (pRNA,snRNA and snoRNA) in the annotation list
  other_family <- c("piRNA", "snRNA", "snoRNA")
  
  if (missing(sncrna_family)) {
    sncrna_family = c(gtsrna_family, mtsrna_family, rsrna_family, ysrna_family)
  } else if(sncrna_family == "tsrna"){
    sncrna_family = c(gtsrna_family, mtsrna_family)
  } else if(sncrna_family == "rsrna"){
    sncrna_family = rsrna_family
  } else if(sncrna_family == "ysrna"){
    sncrna_family = ysrna_family
  } else if (sncrna_family == "mirna"){
    sncrna_family = mirna_family
  } else if (sncrna_family == "other"){
    sncrna_family = other_family
  } else{
    sncrna_family = c(gtsrna_family, mtsrna_family, rsrna_family, ysrna_family, mirna_family, other_family)
  }
  
  #### max_num : number (default 1000) of top species for each sncRNA family to be considered for analysis
  max_num = as.numeric(top_species)
  
  #### row_mean : mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
  row_mean = as.integer(row_mean)
  
  #### p_correction_method : adjustment method ( default : the Bonferroni correction) for correcting or adjusting the p-values
  #p_correction_method = "bonferroni"
  #padj_method = as.character(padj_method)
  if (missing(padj_method)) {
    padj_method = "bonferroni"
  } else if(padj_method == "bonferroni"){ 
    padj_method = "bonferroni"
  } else if(padj_method == "BH"){ 
    padj_method = "BH"
  } else{
    print("Error: use either bonferroni or BH for correcting or adjusting the p-values.")
    return(NULL)
  }
  
  ########## Part2: Function body starts to perform differential expression analysis from processed inputs ##########
  
  if (ncol(e) != length(cl))
  {
    print("Error: column number doesn't match. Please correct the sample numbers.")
    return(NULL)
  }
  
  e = e[rowMeans(e) > row_mean,]   #### only retain the sncRNA species with mean RPM > row_mean
  
  z = match(rownames(e), anno$Sequence)
  if (length(z[is.na(z)]) > 0)
  {
    print("Error: annotation cannot match the input data")
    return(NULL)
  }
  anno = anno[z,]
  
  sncrna_family_species_n = rep(NA, length(sncrna_family))
  
  if (length(sncrna_family) < 1) { 
    print ("There are no species annotated for the chosen sncRNA family in the input matrix at the given threshold")
    return(NULL)
    }
  else {
  for (k in 1:length(sncrna_family))
  {
    if (grepl("rRNA", sncrna_family[k])) anno_k = anno[anno$Annotation==sncrna_family[k],]
    
    else if(grepl("mir", sncrna_family[k], ignore.case = TRUE)) anno_k = anno[grepl(paste0("\\b", sncrna_family[k], "\\b"), anno$Annotation), ]
    else if(grepl("let", sncrna_family[k])) anno_k = anno[grepl(paste0("\\b", sncrna_family[k], "\\b"), anno$Annotation), ]
    
    else anno_k = anno[grepl(sncrna_family[k], anno$Annotation),]
    sncrna_family_species_n[k] = nrow(anno_k)
  }
   }
  sncrna_family
  sncrna_family_species_n
  
  
  sncrna_family = sncrna_family[sncrna_family_species_n >= 2]
  t = rep(NA, length(sncrna_family))
  p = rep(NA, length(sncrna_family))
  
  for (k in 1:length(sncrna_family))
  {
    if (grepl("rRNA", sncrna_family[k])) is.sncrna_family = (anno$Annotation==sncrna_family[k])
    
    else if(grepl("mir", sncrna_family[k], ignore.case = TRUE)) is.sncrna_family = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
    else if(grepl("let", sncrna_family[k])) is.sncrna_family = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
    
    else is.sncrna_family = grepl(sncrna_family[k], anno$Annotation)
    e_tmp = e[is.sncrna_family,]
    
    if (nrow(e_tmp) > max_num)
    {
      e_tmp = e_tmp[order(rowMeans(e_tmp), decreasing=T),]
      e_tmp = e_tmp[1:max_num,]
    }
    
    data = data.frame(rep(rownames(e_tmp), ncol(e_tmp)), unlist(log10(e_tmp+0.1)), rep(cl, rep(nrow(e_tmp), ncol(e_tmp))))
    colnames(data) = c("RNA", "expr", "cl")
    
    lm_summary = summary(lm(data$expr ~ data$RNA + data$cl))
    t[k] = lm_summary$coefficients["data$cl",][3]
    p[k] = lm_summary$coefficients["data$cl",][4]
  }
  adjusted_p = p.adjust(p, method=padj_method)    #### Correcting or adjusting the p-values with padj_method (i.e. the Bonferroni correction (default) or by BH (Benjamini & Hochberg))
  output = data.frame(sncrna_family, t, p, adjusted_p)      #### writing final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
  
  return(output)
}
########## FUSION_ms function ends ##########


#' FUSION_msmc - Differential expression analysis of sncRNA families in multiple samples data with multiple conditions using expression matrix
#'
#' @param a a matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Rest of the columns provides RPM or expression values from different samples under study. Sequence (or ID) must be unique.
#' @param cl a file specifying multiple sample conditions in a comma separated format (such as for 3 different conditions with 6 samples each the file will contains the input: 1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3). Note: one can specify the condition in any order depending on the order of  the samples in the matrix, such as : “2,2,2,1,1,1,3,3,3” or “1,2,1,2,1,2,1,2,1,2” or “1,1,1,2,2,2,1,1,1,3,3,3”.
#' @param row_mean  mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
#' @param sncrna_family list of sncRNA families to be analysis for expression analysis study. Use "tsrna", "rsrna", "ysrna", or "mirna" for tsRNAs, rsRNAs, ysRNAs, or mirna, respectively. Use "other" for 'pRNA,snRNA and snoRNA'. For all, use any letter or number, eg. "a","b", "c", 1, 2, 3. By default (i.e. if no option specified) it will search for tryRNAs (tsRNAs, rsRNAs, and ysRNAs)
#' @param top_species number (default 1000) of top species for each sncRNA family to be considered for analysis. It can help to reduce the run time of the analysis for the families (such as rsrna families) with large number of species. If time is not a concern, use high values such as 5000, 10000, etc.
#' @param padj_method adjustment methods for correcting p-values. Use either "bonferroni" or "BH", where BH stands for Benjamini & Hochberg. By default (i.e. if no option specified) it will run for "bonferroni". 
#' @return It will return a final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
#'
#' @examples
#' # Note:  After installation, one can find the example files in "../FUSION/extdata/"
#' # To run on your own matrix file, provide the full path as : 
#' # FUSION(a = "/path/to/your_matrix.txt", cl = "/path/to/your_condition.txt")
#' example_matrix_cl <- system.file("extdata", "example_matrix_cl.txt", package = "FUSION")
#' example_condition1 <- system.file("extdata", "example_condition1.txt", package = "FUSION")
#' example_condition2 <- system.file("extdata", "example_condition2.txt", package = "FUSION")
#' example_condition3 <- system.file("extdata", "example_condition3.txt", package = "FUSION")
#' FUSION_msmc(a = example_matrix_cl, cl = example_condition1, row_mean = 1, top_species = 5000)   
#' # For running differential expression analysis on example_matrix_cl.txt with 18 samples as per
#' # the conditions specified (i.e., 3 different conditions with 6 samples each) in the file 
#' # example_condition1.txt at row_mean threshold of 1 for 5000 top_species for 
#' # sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
#' FUSION_msmc(a = example_matrix_cl, cl = example_condition2, row_mean = 1, top_species = 5000)   
#' # For running differential expression analysis on example_matrix_cl.txt with 18 samples as per
#' # the conditions specified (i.e., 4 different conditions with condition 1, 2, and 3 having five
#' # samples each, while last three samples are representing the condition 4) in the file 
#' # example_condition2.txt at row_mean threshold of 1 for 5000 top_species for 
#' # sncRNA families (tsRNAs, rsRNAs and ysRNAs) (default);
#' FUSION_msmc(a = example_matrix_cl, cl = example_condition3, row_mean = 10, sncrna_family = "tsrna") 
#' # For running differential expression analysis on example_matrix_cl.txt with 18 samples as per
#' # the conditions specified (i.e., 2 different conditions with condition 1, and 2 having  10 and 
#' # 8 samples, respectively) in the file example_condition3.txt at row_mean threshold of 10 for 
#' # default (1000) top_species for tsRNA families. 
#' # Note:  If you want to save the terminal/console output to a file, use sink() command. 
#' # e.g.,:  options(max.print = 1e6); sink("~/output.txt"); 
#' # FUSION_ms(a = "./extdata/example_matrix1.txt", S1 = 10, S2 = 16); sink()
#' @export

  FUSION_msmc <- function(a, cl, row_mean = getOption("row_mean", "0.1"), sncrna_family = c("tsrna", "rsrna", "ysrna", "mirna", "other"), top_species = getOption("max_num", "1000"), padj_method = c("bonferroni", "BH")){
  #### FUSION_msmc function body starts ####
  args = commandArgs(trailingOnly=TRUE)
  
  
  ########## Part1: Processing of input matrix and different parameters to be used in the function ##########
  #### a: a user-provided matrix with first and second column as "Sequence" (or ID) and "Annotation", respectively. Rest of the columns provides RPM or expression values from different samples under study. Sequence (or ID) must be unique.####
  a <- read.delim(a, header = TRUE)
  a <- as.data.frame(a)
  
  #### e and anno are sub-matrix generated from the user-provided input matrix (a)
  #### e : An expression matrix with first column as Sequence (or ID) and rest of the columns with expression values
  #### anno :  Sequence (or ID) and their Annotation information
  
  anno = a[,1:2]  #### Separated the columns with unique Sequence (or ID) and their Annotation information
  e = a[,3:ncol(a)] #### Separated the columns with expression values
  rownames(e) = a[,1]
  
  #### cl: a vector specifying samples from different conditions (1 and 2) provided as file
  cl <- scan(cl, sep = ",")  
  # Custom message
  cat("Successfully read", length(cl), "samples data from the file.\n")
  
  # View the vector
  cl <- as.numeric(cl)  # Convert the captured output back to a numeric vector
  
  ##### Extracting unique annotation list from input file
  uni_annotation = unique(unlist(strsplit(format(anno[,2], justify="none"), ";")))
  
  #### Extracting genomic tRNAs in the annotation list
  tmp1 = uni_annotation[grepl("mature-tRNA", uni_annotation)]
  gtsrna_family = sort(unique(substr(tmp1, 1, 19)))
  gtsrna_family[gtsrna_family=="mature-tRNA-iMet-CA"] = "mature-tRNA-iMet-CAT"
  
  #### Extracting mitochondrial tRNAs in the annotation list
  tmp2 = uni_annotation[grepl("mature-mt_tRNA", uni_annotation)]
  mtsrna_family = sort(unique(substr(tmp2, 1, 22)))
  
  #### Extracting rRNAs in the annotation list
  #rsrna_family = c("4.5S-rRNA", "5S-rRNA", "5.8S-rRNA", "12S-rRNA", "16S-rRNA", "18S-rRNA", "28S-rRNA", "45S-rRNA")
  tmp3 = uni_annotation[grepl("S-rRNA", uni_annotation)]
  rsrna_family = sort(unique(substr(tmp3, 1, 9)))
  
  #### Extracting YRNAs in the annotation list
  #ysrna_family = c("RNY1-YRNA", "RNY3-YRNA", "RNY4-YRNA", "RNY5-YRNA")
  tmp4 = uni_annotation[grepl("-YRNA", uni_annotation)]
  ysrna_family = sort(unique(substr(tmp4, 1, 9)))
  
  #### Extracting miRNAs in the annotation list
  mirna_patterns <- c("-lin", "-let", "-mir", "-miR", "-MIR")
  tmp5 = uni_annotation[grepl(paste(mirna_patterns, collapse='|'), uni_annotation)]
  mirna_family = sort(unique(substr(tmp5, 1, 16)))
  
  #### Extracting other sncRNAs (pRNA,snRNA and snoRNA) in the annotation list
  other_family <- c("piRNA", "snRNA", "snoRNA")
  
  if (missing(sncrna_family)) {
    sncrna_family = c(gtsrna_family, mtsrna_family, rsrna_family, ysrna_family)
  } else if(sncrna_family == "tsrna"){
    sncrna_family = c(gtsrna_family, mtsrna_family)
  } else if(sncrna_family == "rsrna"){
    sncrna_family = rsrna_family
  } else if(sncrna_family == "ysrna"){
    sncrna_family = ysrna_family
  } else if (sncrna_family == "mirna"){
    sncrna_family = mirna_family
  } else if (sncrna_family == "other"){
    sncrna_family = other_family
  } else{
    sncrna_family = c(gtsrna_family, mtsrna_family, rsrna_family, ysrna_family, mirna_family, other_family)
  }
  
  #### max_num : number (default 1000) of top species for each sncRNA family to be considered for analysis
  max_num = as.numeric(top_species)
  
  #### row_mean : mean RPM (default 0.1) threshold to retain the sncRNA species (rows) in the matrix
  row_mean = as.integer(row_mean)
  
  #### p_correction_method : adjustment method ( default : the Bonferroni correction) for correcting or adjusting the p-values
  #p_correction_method = "bonferroni"
  #padj_method = as.character(padj_method)
  if (missing(padj_method)) {
    padj_method = "bonferroni"
  } else if(padj_method == "bonferroni"){ 
    padj_method = "bonferroni"
  } else if(padj_method == "BH"){ 
    padj_method = "BH"
  } else{
    print("Error: use either bonferroni or BH for correcting or adjusting the p-values.")
    return(NULL)
  }
  
  ########## Part2: Function body starts to perform differential expression analysis from processed inputs ##########
  
  if (ncol(e) != length(cl))
  {
    print("Error: column number doesn't match. Please correct the sample numbers.")
    return(NULL)
  }
  
  e = e[rowMeans(e) > row_mean,]   #### only retain the sncRNA species with mean RPM > row_mean
  
  z = match(rownames(e), anno$Sequence)
  if (length(z[is.na(z)]) > 0)
  {
    print("Error: annotation cannot match the input data")
    return(NULL)
  }
  anno = anno[z,]
  
  sncrna_family_species_n = rep(NA, length(sncrna_family))
  
  if (length(sncrna_family) < 1) { 
    print ("There are no species annotated for the chosen sncRNA family in the input matrix at the given threshold")
    return(NULL)
  }
  else {
    
  for (k in 1:length(sncrna_family))
  {
    if (grepl("rRNA", sncrna_family[k])) anno_k = anno[anno$Annotation==sncrna_family[k],]
    
    else if(grepl("mir", sncrna_family[k], ignore.case = TRUE)) anno_k = anno[grepl(paste0("\\b", sncrna_family[k], "\\b"), anno$Annotation), ]
    else if(grepl("let", sncrna_family[k])) anno_k = anno[grepl(paste0("\\b", sncrna_family[k], "\\b"), anno$Annotation), ]
    
    else anno_k = anno[grepl(sncrna_family[k], anno$Annotation),]
    sncrna_family_species_n[k] = nrow(anno_k)
  }
  }
  sncrna_family
  sncrna_family_species_n
  
  
  sncrna_family = sncrna_family[sncrna_family_species_n >= 2]
  t = rep(NA, length(sncrna_family))
  p = rep(NA, length(sncrna_family))
  
  for (k in 1:length(sncrna_family))
  {
    if (grepl("rRNA", sncrna_family[k])) is.sncrna_family = (anno$Annotation==sncrna_family[k])
    
    else if(grepl("mir", sncrna_family[k], ignore.case = TRUE)) is.sncrna_family = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
    else if(grepl("let", sncrna_family[k])) is.sncrna_family = grepl((paste0("\\b", sncrna_family[k], "\\b")), anno$Annotation)
    
    else is.sncrna_family = grepl(sncrna_family[k], anno$Annotation)
    e_tmp = e[is.sncrna_family,]
    
    if (nrow(e_tmp) > max_num)
    {
      e_tmp = e_tmp[order(rowMeans(e_tmp), decreasing=T),]
      e_tmp = e_tmp[1:max_num,]
    }
    
    data = data.frame(rep(rownames(e_tmp), ncol(e_tmp)), unlist(log10(e_tmp+0.1)), rep(cl, rep(nrow(e_tmp), ncol(e_tmp))))
    colnames(data) = c("RNA", "expr", "cl")
    
    lm_summary = summary(lm(data$expr ~ data$RNA + data$cl))
    t[k] = lm_summary$coefficients["data$cl",][3]
    p[k] = lm_summary$coefficients["data$cl",][4]
  }
  adjusted_p = p.adjust(p, method=padj_method)    #### Correcting or adjusting the p-values with padj_method (i.e. the Bonferroni correction (default) or by BH (Benjamini & Hochberg))
  output = data.frame(sncrna_family, t, p, adjusted_p)      #### writing final output in a data-frame with t-statistics, P_value, and adjusted P_value for each sncRNA family chosen for analysis
  
  return(output)
}
########## FUSION_msmc function ends ##########
