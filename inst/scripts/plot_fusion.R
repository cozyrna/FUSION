library(Biostrings)
library(seqinr)

#input_rpm_file: RPM data of the paired samples
#parental_rna_file: fasta file of the parental RNA
#top_rna_num: number of RNA species to be visualized (default 100)
#pos_coord: start and end positions within the parental RNA to be visualized
#max.mismatch: maximum number of mismatch (default 1)

# === Set input files ===
# By default, use example files bundled with the package:
parental_rna_file <- system.file("extdata", "human_rRNA_28S.fa", package = "FUSION")
input_rpm_file <- system.file("extdata", "example_visualization_data.txt", package = "FUSION")

# === To use your own data, uncomment and edit the lines below ===
# parental_rna_file <- "/full/path/to/your/input.fa"
# input_rpm_file <- "/full/path/to/your/input_rpm.txt"

# === Load parental RNA sequence ===
parental_rna <- unlist(read.fasta(parental_rna_file, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE))

# === Plotting function ===
plot.fusion_ps <- function(input_rpm_file, parental_rna, top_rna_num = 100, pos_coord = NULL, max.mismatch = 1)
{
  a = read.delim(input_rpm_file, comment.char="#")
  e = a[,2:3]
  rownames(e) = a$Sequence
  e = e[order(rowMeans(e), decreasing=T),][1:top_rna_num,]
  
  ref_rna = Biostrings::DNAStringSet(parental_rna)
  start = rep(NA, nrow(e))
  len = rep(NA, nrow(e))
  for (i in 1:nrow(e))
  {
    input_frag = rownames(e)[i]
    hit = as.data.frame(unlist(vmatchPattern(pattern = input_frag, subject = ref_rna, max.mismatch = max.mismatch)))
    if (nrow(hit) < 1) next
    start[i] = hit[1,1]
    len[i] = hit[1,3]
  }
  
  if (is.null(pos_coord)) xlim = c(0, nchar(parental_rna))
  else xlim = pos_coord
  
  plot(1, xlim=xlim, ylim=range(e), xlab="Position (nt)", ylab="RPM", main="", type="n", axes=F, log="y")
  for (i in 1:nrow(e))
  {
    lines(c(start[i], start[i]+len[i]-1), c(e[i,1], e[i,1]), col="slateblue", lwd=1)
    lines(c(start[i], start[i]+len[i]-1), c(e[i,2], e[i,2]), col="red", lwd=1)
    lines(c(start[i], start[i]), c(e[i,1], e[i,2]), col="grey60", lwd=0.5, lty=2)
    lines(c(start[i]+len[i]-1, start[i]+len[i]-1), c(e[i,1], e[i,2]), col="grey60", lwd=0.5, lty=2)
    points(c(start[i], start[i]+len[i]-1), c(e[i,1], e[i,1]), pch=19, col="slateblue", cex=0.5)
    points(c(start[i], start[i]+len[i]-1), c(e[i,2], e[i,2]), pch=19, col="red", cex=0.5)
  }
  axis(1)
  axis(2)
}


# === Generate the plots ===
pdf("Parent_RNA_plot.pdf")
#plot the whole parental RNA
plot.fusion_ps(input_rpm_file, parental_rna)
#plot the zoomed-in view of the specified region within the parental RNA
plot.fusion_ps(input_rpm_file, parental_rna, pos_coord = c(640, 880)) #mention the coordinates here for zoomed-in positions
dev.off()
cat("File 'Parent_RNA_plot' is written to:", getwd(), "\n\n")

