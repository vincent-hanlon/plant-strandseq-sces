
# Example breakpointR run for examining Strand-seq libraries aligned to a reference genome.
# Chromosome names may differ
# Run this in a directory of single-cell BAM files for Strand-seq libraries


library(breakpointR)

chr <- c("CM200447.1","CM200448.1","CM200449.1","CM200450.1","CM200451.1","CM200452.1","CM200453.1","CM200454.1","CM200455.1","CM200456.1","CM200457.1","CM200458.1","CM200459.1","CM200460.1","CM200461.1","CM200462.1","CM200463.1")

pdf()

breakpointr(inputfolder="./", outputfolder="./BPR_output", pairedEndReads=TRUE,numCPU=24,windowsize=1000000,binMethod="size", 
	background=0.15, chromosomes=chr, maskRegions=NULL)

dev.off()
