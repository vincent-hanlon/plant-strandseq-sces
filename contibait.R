# R code for the pseudocontig analysis for Moxley et al.
# https://www.bioconductor.org/packages/release/bioc/html/contiBAIT.html
# Edited from the contiBAIT tutorial/vignette for PE sunflower libraries
# I just source()'d this script from R in a directory containing only the good-quality sunflower libraries
# Use good-quality libraries only!

library(contiBAIT)
library(parallel)


# Splitting sunflower genome into 1 Mb pseudocontigs for re-scaffolding demo
bamFileList <- list.files(path = ".", pattern=".bam$")
dividedChr <- makeChrTable(bamFileList[1], splitBy=1000000)

# Getting strand states of 1 Mb bins; excluding bins with less than 15 reads
# (rather than the default, 10)
strandFrequencyList <- strandSeqFreqTable(bamFileList, 
	filter=dividedChr,
	qual=10, 
	pairedEnd=TRUE,
	BAITtables=TRUE)

strandFreq <- strandFrequencyList[[1]]
readCounts <- strandFrequencyList[[2]]
strandFreq[which(readCounts < 15)] <- NA 

pdf("strandseq_quality.pdf")
plotWCdistribution(strandFreq)
dev.off()

# Building the strandStateMatrix
strandStateMatrix <- preprocessStrandTable(
	strandFreq, 
	lowQualThreshold=0.8)

WCMatrix <- strandStateMatrix[[1]]

# Clustering contigs into chromosomes and reorienting backwards ones
clusteredContigs <- clusterContigs(WCMatrix, randomise=TRUE, recluster=12, randomSeed=32459)
reorientedMatrix <- reorientAndMergeLGs(clusteredContigs, WCMatrix, cluster=12)

pdf("linkage_group_distances.pdf")
plotLGDistances(reorientedMatrix[[3]], reorientedMatrix[[1]])
dev.off()

pdf("contig_order_heatmap.pdf")
# Ordering the contigs within chromosomes
# Some debate with the authors of this package whether the reoriented matrix should be provided to orderAllLinkageGroups, or whether the original WCMatrix should be
# It looks a little better when the reoriented matrix is used, which maybe makes sense, but that's not what's in the vignette.
contigOrder <- orderAllLinkageGroups(reorientedMatrix[[3]],
	reorientedMatrix[[1]],
	strandFreq,
	readCounts,
	saveOrdered=TRUE, nProcesses=12, randomAttempts=100)
dev.off()

# More plotting for QC
pdf("contig_order_scatter.pdf")
plotContigOrder(contigOrder[[1]])
dev.off()

WatsonFreqList <- strandFrequencyList[[3]]
CrickFreqList <- strandFrequencyList[[4]]

pdf("plot_chrs_all_libs.pdf")
ideogramPlot(WatsonFreqList,
	CrickFreqList,
	dividedChr,
	orderFrame=contigOrder[[1]],
	plotBy='chr',
	orientationData = reorientedMatrix[[2]])
dev.off()
write.table(as.data.frame(contigOrder[[1]]), "contig_order.txt", sep="\t", quote=F, col.names=F,row.names=F)
r <- reorientedMatrix[[2]];write.table(data.frame(seqnames(r), start(r), end(r), strand(r)), "contig_orientation.txt", quote=F, sep="\t", col.names=F, row.names=F)
