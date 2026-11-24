# Ran this on all good-quality Strand-seq libraries (single-cell BAM files).
# The output file breakPointSummary.txt corresponds to sces_unprocessed.txt
# We then removed strand switches (lines) that were (i) in borderline-quality libraries listed in spurious_high_SCE_libraries.txt (ii) had ww-cc or cc-ww strand-switches (iii) were within 5 Mb of each other within the same library (midpoint to midpoint)

library(breakpointR)

	chr <- c("Ha412HOv3Chr01","Ha412HOv3Chr02","Ha412HOv3Chr03","Ha412HOv3Chr04","Ha412HOv3Chr05","Ha412HOv3Chr06","Ha412HOv3Chr07","Ha412HOv3Chr08","Ha412HOv3Chr09","Ha412HOv3Chr10","Ha412HOv3Chr11","Ha412HOv3Chr12","Ha412HOv3Chr13","Ha412HOv3Chr14","Ha412HOv3Chr15","Ha412HOv3Chr16","Ha412HOv3Chr17")
	mask <- 'bpr_residual_strand_switch_mask.bed'


pdf()

breakpointr(inputfolder="./", outputfolder="./BPR_output_5MB_masked", pairedEndReads=TRUE,numCPU=24,windowsize=5000000,binMethod="size", 
	background=0.15, chromosomes=chr, maskRegions=mask)

dev.off()
