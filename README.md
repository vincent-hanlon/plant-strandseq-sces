This repository contains scripts relating to the article

 > **Plant Strand-seq resolves fine-scale genome structure and reveals genome-wide patterns of DNA damage** \
 > Dylan R. Moxley*, Vincent C.T. Hanlon*, Cassandra Elphinstone, Eric Gerardo González Segovia, Yanni Wang, Daniel Chan, Loren H. Rieseberg, Peter M. Lansdorp, Marco Todesco. (* equal contributions) \

### Sister chromatid exchange (SCE) analysis
Running [breakpointr_sce_calling.R](SCEs/breakpointr_sce_calling.R) on a set of single-cell BAM files from good-quality Strand-seq libraries aligned to the HA412HOv3 sunflower genome assembly produces a callset of strand-switches which can be further filtered to generate SCE calls. [breakpointr_sce_calling.R](SCEs/breakpointr_sce_calling.R) uses [breakpointR](https://github.com/daewoooo/breakpointR) with large 5 Mb bins (to identify "big" events that are likely SCEs) and masks a set of recurrent strand-switches that are not SCEs in [bpr_residual_strand_switch_mask.bed](SCEs/bpr_residual_strand_switch_mask.bed).

We processed the resulting callset (output as 'breakPointSummary.txt' but renamed in this repo to [sces_unprocessed.txt](SCEs/sces_unprocessed.txt)) by removing: 
1. strand switches in the borderline-poor-quality libraries listed in [spurious_high_SCE_libraries.txt](SCEs/spurious_high_SCE_libraries.txt).
2. "ww-cc" and "cc-ww" strand switches (as SCEs should affect just one homolog).
3. events within 5 Mb of each other within the same cell and on the same chromosome (measured midpoint to midpoint).

The final SCE callset is in [sces_processed.txt](SCEs/sces_processed.txt). A BED file with regions in which SCEs cannot be reliably called, because of poor read mapping at telomeres and because there must be half of a 5 Mb bin on either side of an SCE to call it, is provided [here](SCEs/sce_chromosome_end_mask.bed)

### Scaffolding demonstration
Again using a set of good-quality Strand-seq libraries, but this time aligned to the HA412HOv2 sunflower genome assembly, we used [contiBAIT](https://doi.org/10.1093/bioinformatics/btx281) to break the assembly into 1 Mb pseudocontigs and re-scaffold them with Strand-seq alone. In R, `source()` the script [contibait.R](contibait.R) from a directory with (only) the good-quality BAMs.

### Misorient analysis
Strand-seq reads encode the orientation of DNA template strands at the previous cell division: so both "backwards" region in a genome assembly and polymorphic inversions generate a region in which Strand-seq reads map in the opposite orientation reads in neighbouring regions. We therefore aligned Strand-seq reads to HA412HOv2 to identify misoriented segmnents within the genome assembly. [InvertypeR](https://github.com/vincent-hanlon/InvertypeR) is run as shown in (invertyper.R)[invertyper.R], pretending the sunflower reads come from a single haploid individual because it is highly homozygous, and because we are interested in assembly misorients that look like inversions in a haploid individual. We must mask some [very large misorients on chromosome 16](invertyper_soft_mask.bed) because they interfere with composite file creation. Misorient calls can be found [here](all_inversions_inc_contiBAIT.txt) with the addition of one very large misorient that InvertypeR did not find but contiBAIT did.

### Miscellaneous
The script [breakpointr_example.R](breakpointr_example.R) gives an example of a [breakpointR](https://github.com/daewoooo/breakpointR) run to examine Strand-seq libraries visually (e.g. for assessing quality, or for checking SCEs).
[This file](good_quality_strandseq_libraries.list) contains a list of Strand-seq libraries that passed manual quality curation ("good quality").
