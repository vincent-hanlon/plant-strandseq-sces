# This is the command that can be used to call misorients from the sunflower libraries using InvertypeR
# More details at https://github.com/vincent-hanlon/InvertypeR
# Note that because this individual is highly homozygous and we don't have SNVs, we're pretending it's haploid so that InvertypeR doesn't try to create an SNV-phased composite file
# Good quality libraries only!

# Generating composite files may not be quick, so it is parallelized 24-fold.
# Note that the soft mask is needed so that the composite files are created correctly despite the large misorients on chr16 (CM200462.1) in the v2 assembly

# Finally, the chromosome names for any particular version of the assembly may need to be altered here and in the soft mask file.

library(invertyper)

chr <- c("CM200447.1","CM200448.1","CM200449.1","CM200450.1","CM200451.1","CM200452.1","CM200453.1","CM200454.1","CM200455.1","CM200456.1","CM200457.1","CM200458.1","CM200459.1","CM200460.1","CM200461.1","CM200462.1","CM200463.1")

i <- invertyper_pipeline(regions_to_genotype=NULL, hard_mask=NULL, soft_mask="invertyper_soft_mask.bed",
    chromosomes=chr, haploid_chromosomes=chr, vcf=NULL, numCPU=24, adjust_method='all', 
    save_composite_files=T, write_browser_files=T, discover_breakpointr_inversions=T)
