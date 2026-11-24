This repository contains scripts relating to the article

 > **Plant Strand-seq resolves fine-scale genome structure and reveals genome-wide patterns of DNA damage** \
 > Dylan R. Moxley*, Vincent C.T. Hanlon*, Cassandra Elphinstone, Eric Gerardo González Segovia, Yanni Wang, Daniel Chan, Loren H. Rieseberg, Peter M. Lansdorp, Marco Todesco. \
 > (* equal contributions)

### Sister chromatid exchange (SCE) analysis
Running [SCEs/breakpointr_sce_calling.R](breakpointr_sce_calling.R) on a set of single-cell BAM files from good-quality Strand-seq libraries aligned to the HA412HOv3 sunflower genome assembly produces a callset of strand-switches which can be further filtered to generate SCE calls. [SCEs/breakpointr_sce_calling.R](breakpointr_sce_calling.R) uses large 5 Mb bins (to identify "big" events that are likely SCEs) and masks a set of recurrent strand-switches that are not SCEs in [SCEs/bpr_residual_strand_switch_mask.bed](bpr_residual_strand_switch_mask.bed).

We processed the resulting callset (output as 'breakPointSummary.txt' but renamed in this repo to [SCEs/sces_unprocessed.txt](sces_unprocessed.txt)) by: 
    (i) Removing strand-switches in the borderline-poor-quality libraries listed in [SCEs/spurious_high_SCE_libraries.txt](spurious_high_SCE_libraries.txt)
    (ii) Removing "ww-cc" and "cc-ww" strand switches (as SCEs should affect just one homolog)
    (iii) Removing events within 5 Mb of each other within the same cell and on the same chromosome (measured midpoint to midpoint).

The final SCE callset is in [SCEs/sces_processed.txt](sces_processed.txt). A BED file with regions in which SCEs cannot be reliably called, because of poor read mapping at telomeres and because there must be half of a 5 Mb bin on either side of an SCE to call it, is provided [SCEs/sce_chromosome_end_mask.bed](here)

