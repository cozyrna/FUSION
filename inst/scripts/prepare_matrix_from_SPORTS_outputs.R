#!/usr/bin/Rscript

#### USage: Rscript prepare_matrix_from_SPORTS_outputs.R file_list_document out_prefix
## 'file_list_document' must contains the full path and name of X_output.txt files from the SPORTS output.
## The order in which the files are listed will determine the sample order in the matrix file. Please ensure the files are arranged accordingly.

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library("tidyverse")

####    merge SPORTS output files to show annotation with read counts across samples - START
# Read all files names in the list
files.Name <- args[1]
conn <- file(files.Name, open="r")
all_lines <- readLines(conn)
lines <- all_lines[which(all_lines != "")]
close(conn)

# Extract the base filenames (without path and without .extension)
base_filenames <- tools::file_path_sans_ext(basename(lines))
  
# Remove the text "_output" from the filenames
base_filenames <- gsub("_output", "", base_filenames)
  
# add other column headers that will precede in the final matrix
final_row <- c("Sequence",	"Length",	"Match_Genome",	"Annotation", base_filenames)

# change name of "Reads" column to keep it unique for each sample
len <-length(lines) 
file_list <- list()
for (i in 1:len){
# Loop through each file
  file_list[[i]] <- read.delim(lines[i], header=TRUE)
  colnames(file_list[[i]])[4] = paste0("Sample_", i)
}
gc()

## merge all data frames into one, only differ by read counts  
oneDF <- file_list %>% reduce(full_join, by = c("Sequence",	"Length",	"Match_Genome",	"Annotation"))
# remove IDs columns
oneDF1 <- oneDF %>% select(-starts_with("ID"))
# replace 'NA" with "0"
oneDF1[is.na(oneDF1)] <- 0
# Rearrange columns
oneDF1 <- oneDF1 %>% relocate(Sample_1, .after = Annotation)
colnames(oneDF1) <- final_row[1:ncol(oneDF1)]

# Convert to data.table
oneDT <- as.data.table(oneDF1)
gc()
# Add Annotation column to keep one annotation per sequence
oneDT[, Annotation := Annotation[1], by = .(Sequence, Length, Match_Genome)]

# Summarize numeric columns (sum them) grouped by Sequence, Length, Match_Genome, and Annotation
oneDF2 <- oneDT[, lapply(.SD, sum), by = .(Sequence, Match_Genome, Annotation), .SDcols = is.numeric]
gc()

# Retain Length for each group (keeping the first value)
oneDF2[, Length := oneDT[, first(Length), by = .(Sequence, Match_Genome)]$V1]

# Reorder columns to match the desired order
setcolorder(oneDF2, c("Sequence", "Length", "Match_Genome", "Annotation", setdiff(names(oneDF2), c("Sequence", "Length", "Match_Genome", "Annotation"))))
gc()
####  merge SPORTS output files to show annotation with read counts across samples - END

#### generate read_count matrix from merged output for downstream analysis - FUSION, DESeq2, etc. - START
# only keep the Sequence, Annotation and Read_counts columns
setDT(oneDF2)
Read_Matrix_DF <- oneDF2[, !c("Length", "Match_Genome"), with = FALSE]
gc()

# write read counts table to a file
data.table::fwrite(Read_Matrix_DF, file= paste0(args[2], "_count-matrix.txt"), sep='\t', row.names=FALSE, quote=FALSE)
gc()
#### generate read_count matrix from merged output for downstream analysis - FUSION, DESeq2, etc.- END

#### generate RPM matrix from count matrix - START
# Consider only the columns with read counts 
Count_data <- Read_Matrix_DF[,-(1:2)]

# Generate Per Million scaling factor and then RPM count per column
scaling_factor <- apply(Count_data, 2, sum)/1000000
RPM_Matrix <- mapply('/', Count_data, scaling_factor)
gc()

# add "Sequence" column to generate RPM counts matrix
RPM_Matrix_DF <- cbind(Read_Matrix_DF[,1], Read_Matrix_DF[,2], RPM_Matrix)

# write RPM counts matrix table to a file
data.table::fwrite(RPM_Matrix_DF, file= paste0(args[2], "_RPM-matrix.txt"), sep='\t', row.names=FALSE, quote=FALSE)
#### generate RPM matrix from count matrix - END
