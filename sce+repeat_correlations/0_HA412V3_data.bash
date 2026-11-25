###################################
# Datasets for repeat, SCE and A/B compartment correlations
# November 2026
####################################

# Datasets
cd /scratch/celphin/repeats/sunflower2/SCE_repeats_corr/

# masking files showing the extent of each telomere to exclude

more telomere_mask_for_5Mb_bins.bed
# Ha412HOv3Chr01  0       2500000
# Ha412HOv3Chr01  147277849       149777849
# Ha412HOv3Chr02  0       4500000
# Ha412HOv3Chr02  173904745       176404745
# Ha412HOv3Chr03  0       4500000
# Ha412HOv3Chr03  171114061       173614061

#-------------------
# repeat annotation
more TEsorter_annotation.ALLchrs.txt
# Name    Chr     Start   End     Repeatclass     Order   Superfamily     Clade   Complete
# TE_00002742_INT Ha412HOv3Chr01  1       859     Copia_LTR_retrotransposon       .       .       .       .
# TE_00003645     Ha412HOv3Chr01  20283   21687   L1_LINE_retrotransposon LINE    unknown unknown unknown
# TE_00003039_LTR Ha412HOv3Chr01  178553  180085  Gypsy_LTR_retrotransposon       .       .       .       .
# TE_00002491_LTR Ha412HOv3Chr01  180086  180211  LTR_retrotransposon     LTR     Gypsy   Tekay   no
# TE_00003438_INT Ha412HOv3Chr01  180219  183660  Gypsy_LTR_retrotransposon       LTR     Gypsy   Tekay   no
# TE_00002755_INT Ha412HOv3Chr01  182880  185003  LTR_retrotransposon     LTR     Gypsy   mixture no
# TE_00002491_LTR Ha412HOv3Chr01  184841  185007  LTR_retrotransposon     LTR     Gypsy   Tekay   no
# TE_00003382_LTR Ha412HOv3Chr01  185008  185654  LTR_retrotransposon     .       .       .       .
# TE_00002491_LTR Ha412HOv3Chr01  187607  187733  LTR_retrotransposon     LTR     Gypsy   Tekay   no
# TE_00003438_INT Ha412HOv3Chr01  187741  191181  Gypsy_LTR_retrotransposon       LTR     Gypsy   Tekay   no
# TE_00002755_INT Ha412HOv3Chr01  190401  192524  LTR_retrotransposon     LTR     Gypsy   mixture no
# TE_00002491_LTR Ha412HOv3Chr01  192362  192528  LTR_retrotransposon     LTR     Gypsy   Tekay   no

# get unique repeat classes
awk 'NR > 1 {print $5}' TEsorter_annotation.ALLchrs.txt | sort -u
# Bel_Pao_LTR_retrotransposon
# CACTA_TIR_transposon
# Copia_LTR_retrotransposon
# CRE_LINE_retrotransposon
# DNA_transposon
# Endogenous_Retrovirus_LTR_retrotransposon
# Gypsy_LTR_retrotransposon
# hAT_TIR_transposon
# helitron
# L1_LINE_retrotransposon
# long_terminal_repeat
# LTR_retrotransposon
# Mutator_TIR_transposon
# pararetrovirus
# PIF_Harbinger_TIR_transposon
# PILE_TIR_transposon
# POLE_TIR_transposon
# repeat_fragment
# repeat_region
# RTE_LINE_retrotransposon
# snRNA
# target_site_duplication
# Tc1_Mariner_TIR_transposon
# tRNA_SINE_retrotransposon

# make new column grouping LTRs and transposons
awk 'NR == 1 {
       print $0 "\tCategory"
       next
     }
     {
       repeat_class = $5
       if (repeat_class == "long_terminal_repeat")
         category = "long_terminal_repeat"
       else if (repeat_class ~ /LTR/)
         category = "LTR"
       else if (repeat_class ~ /TIR/)
         category = "TIR"
       else
         category = "OTHER"

       print $0 "\t" category
     }' TEsorter_annotation.ALLchrs.txt > TEsorter_annotation.ALLchrs_categorized.txt

more TEsorter_annotation.ALLchrs_categorized.txt
# Name    Chr     Start   End     Repeatclass     Order   Superfamily     Clade   Complete        Category
# TE_00002742_INT Ha412HOv3Chr01  1       859     Copia_LTR_retrotransposon       .       .       .       .       LTR
# TE_00003645     Ha412HOv3Chr01  20283   21687   L1_LINE_retrotransposon LINE    unknown unknown unknown OTHER
# TE_00002742_INT Ha412HOv3Chr01  24982   25763   Copia_LTR_retrotransposon       .       .       .       .       LTR
# TE_00002210_LTR Ha412HOv3Chr01  34875   34960   LTR_retrotransposon     .       .       .       .       LTR
# TE_00002651     Ha412HOv3Chr01  47542   57307   repeat_region   .       .       .       .       OTHER
# TE_00003087_LTR Ha412HOv3Chr01  145819  145995  LTR_retrotransposon     .       .       .       .       LTR
# TE_00002363_LTR Ha412HOv3Chr01  146304  147353  LTR_retrotransposon     .       .       .       .       LTR
# TE_00003746     Ha412HOv3Chr01  147354  151147  repeat_fragment .       .       .       .       OTHER
# TE_00003341     Ha412HOv3Chr01  151143  165073  repeat_region   .       .       .       .       OTHER
# TE_00002491_LTR Ha412HOv3Chr01  180086  180211  LTR_retrotransposon     LTR     Gypsy   Tekay   no      LTR
# TE_00003438_INT Ha412HOv3Chr01  180219  183660  Gypsy_LTR_retrotransposon       LTR     Gypsy   Tekay   no      LTR
# TE_00002755_INT Ha412HOv3Chr01  182880  185003  LTR_retrotransposon     LTR     Gypsy   mixture no      LTR
# TE_00002491_LTR Ha412HOv3Chr01  184841  185007  LTR_retrotransposon     LTR     Gypsy   Tekay   no      LTR
# TE_00003382_LTR Ha412HOv3Chr01  185008  185654  LTR_retrotransposon     .       .       .       .       LTR

#-------------------
# SCEs
more breakPointSummary_switchesMin5Mb_verygoodqualonly.txt
# seqnames start end CI.start CI.end genoT filenames
# Ha412HOv3Chr02 107213485 107328477 106838610 107383419 wc-ww sunflower-22sep29-72h-r03-c40_S4.processed.bam.RData.1
# Ha412HOv3Chr04 146594784 146787846 146275721 146909267 wc-cc sunflower-22sep29-72h-r03-c40_S4.processed.bam.RData.2
# Ha412HOv3Chr07 141451599 141463116 141336065 141713972 wc-ww sunflower-22sep29-72h-r03-c40_S4.processed.bam.RData.3
# Ha412HOv3Chr10 45034752 45047738 44946579 45139989 wc-ww sunflower-22sep29-72h-r03-c40_S4.processed.bam.RData.4

#-------------------
#A/B compartments 0.5Mbp
more Ha412_ABcompval_0.5Mbp.txt
# Ha412HOv3Chr01  1       500000  A       0.01978703954483692     .
# Ha412HOv3Chr01  500001  1000000 B       -0.041127998694391224   .
# Ha412HOv3Chr01  1000001 1500000 B       -0.06192246843243125    .
# Ha412HOv3Chr01  1500001 2000000 B       -0.04624612894858877    .
# Ha412HOv3Chr01  2000001 2500000 B       -0.060601376388977      .
# Ha412HOv3Chr01  2500001 3000000 B       -0.07557090533226578    .
# Ha412HOv3Chr01  3000001 3500000 B       -0.012511783968507743   .

#---------------------
# gene annotations
more Ha412H0V3_genes.bed
# Ha412HOv3Chr02	131799236	131801282	g6195
# Ha412HOv3Chr02	131801317	131802048	g6196
# Ha412HOv3Chr02	132042412	132043791	g6197
# Ha412HOv3Chr02	132063575	132065173	g6198
# Ha412HOv3Chr02	132116416	132117201	g6199
# Ha412HOv3Chr02	132129631	132129837	g6200
# Ha412HOv3Chr02	132132317	132134783	g6201

#--------------------------
# methylation
more HA412HO.5mC.windows_100000.circos.bed
# Ha412HOv3Chr01  0       100000  93.98597
# Ha412HOv3Chr01  100000  200000  96.3578
# Ha412HOv3Chr01  200000  300000  94.24204
# Ha412HOv3Chr01  300000  400000  90.39728
# Ha412HOv3Chr01  400000  500000  93.661545
# Ha412HOv3Chr01  500000  600000  93.771835
# Ha412HOv3Chr01  600000  700000  92.097855
# Ha412HOv3Chr01  700000  800000  92.65173
# Ha412HOv3Chr01  800000  900000  88.89797
# Ha412HOv3Chr01  900000  1000000 91.67708
# Ha412HOv3Chr01  1000000 1100000 87.08699
# Ha412HOv3Chr01  1100000 1200000 92.4267
# Ha412HOv3Chr01  1200000 1300000 91.11026
# Ha412HOv3Chr01  1300000 1400000 83.52785
# Ha412HOv3Chr01  1400000 1500000 93.798294

#--------------------------
# short read mapping counts
more strandseq_readcounts_in_1mb_bins_v3.bed
# Ha412HOv3Chr01  0       1000000 15128
# Ha412HOv3Chr01  1000000 2000000 20599
# Ha412HOv3Chr01  2000000 3000000 25493
# Ha412HOv3Chr01  3000000 4000000 14702
# Ha412HOv3Chr01  4000000 5000000 14370
# Ha412HOv3Chr01  5000000 6000000 19730
# Ha412HOv3Chr01  6000000 7000000 22861
# Ha412HOv3Chr01  7000000 8000000 21101
# Ha412HOv3Chr01  8000000 9000000 21816
# Ha412HOv3Chr01  9000000 10000000        21344
# Ha412HOv3Chr01  10000000        11000000        17432
# Ha412HOv3Chr01  11000000        12000000        14835
# Ha412HOv3Chr01  12000000        13000000        14964
# Ha412HOv3Chr01  13000000        14000000        10048
# Ha412HOv3Chr01  14000000        15000000        9264
# Ha412HOv3Chr01  15000000        16000000        6922

##################################################
