library(Biostrings)
library(seqinr)

#input_rpm_file: RPM data of the paired samples
#parental_rna_file: fasta file of the parental RNA
#top_rna_num: number of RNA species to be visualized (default 100)
#pos_coord: start and end positions within the parental RNA to be visualized
#max.mismatch: maximum number of mismatch (default 1)

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

pdf(paste0("./Parent_RNA_plot.pdf"))  ##add desire output file name with path
parental_rna = unlist(read.fasta("./extdata/human_rRNA_28S.fa", seqtype="DNA", as.string=T, forceDNAtolower=F))   ##add the proper path of the file "human_rRNA_28S.fa"  
#plot the whole parental RNA
plot.fusion_ps("./extdata/example_visualization_data.txt", parental_rna)    ##add the proper path of the file "example_visualization_data.txt"
#plot the zoomed-in view of the specified region within the parental RNA
plot.fusion_ps("./extdata/example_visualization_data.txt", parental_rna, pos_coord=c(640, 880))      ##add the proper path of the file "example_visualization_data.txt"
dev.off()
cat("File 'Parent_RNA_plot' is written to:", getwd(), "\n\n")
