###################################
# Plots of repeat, SCE and A/B compartment along each chromosome
# November 2026
####################################

# rorqual1
tmux new-session -s sunflower
tmux attach-session -t sunflower

cd /scratch/celphin/repeats/sunflower2/SCE_repeats_corr/

module load StdEnv/2020
module load r/4.3.1
R

library(glmmTMB)
library(lme4)
library(lmerTest)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

win_size <- 2000000
fout <- "master_final.csv"

master <- read.csv(fout, header = TRUE)

# check column names
names(master)

###################################
 #Plot boxplot of significant factors for each SCEs and AB

master$sce_present <- as.factor(master$sce_present)
master$sce_count <- as.factor(master$sce_count)
master$compartment <- as.factor(master$compartment)

# Genes
boxplot_sce_ltr <- ggplot(master, aes(x = sce_present, y = frac_gene)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  theme_classic() 
ggsave("frac_gene_by_sce_presence_boxplot_2Mbp_5Mbpmask.png", boxplot_sce_ltr, width = 5, height = 4, dpi = 300)
ggsave("frac_gene_by_sce_presence_boxplot_2Mbp_5Mbpmask.pdf", boxplot_sce_ltr)
# warnings just about masked windows

# A/B compartments
boxplot_sce_ltr <- ggplot(master, aes(x = sce_present, y = comp_val)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  theme_classic() 
ggsave("compartment_by_sce_presence_boxplot_2Mbp_5Mbpmask.png", boxplot_sce_ltr, width = 5, height = 4, dpi = 300)
ggsave("compartment_by_sce_presence_boxplot_2Mbp_5Mbpmask.pdf", boxplot_sce_ltr)


# LTR
boxplot_sce_ltr <- ggplot(master, aes(x = sce_present, y = frac_LTR)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  theme_classic() 
ggsave("LTR_by_sce_presence_boxplot_2Mbp_5Mbpmask.png", boxplot_sce_ltr, width = 5, height = 4, dpi = 300)
ggsave("LTR_by_sce_presence_boxplot_2Mbp_5Mbpmask.pdf", boxplot_sce_ltr)


#---------------
# SCE count/AB compartment
boxplot_sce_ltr <- ggplot(master, aes(x = sce_count, y = comp_val)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  theme_classic() 
ggsave("compartment_by_sce_count_boxplot_2Mbp_5Mbpmask.png", boxplot_sce_ltr, width = 5, height = 4, dpi = 300)
ggsave("compartment_by_sce_count_boxplot_2Mbp_5Mbpmask.pdf", boxplot_sce_ltr)


# SCE count/LTR
boxplot_sce_ltr <- ggplot(master, aes(x = sce_count, y = frac_LTR)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  theme_classic() 
ggsave("frac_LTR_by_sce_count_boxplot_2Mbp_5Mbpmask.png", boxplot_sce_ltr, width = 5, height = 4, dpi = 300)
ggsave("frac_LTR_by_sce_count_boxplot_2Mbp_5Mbpmask.pdf", boxplot_sce_ltr)


######################################
# Plot of frac_LTR for each chromosome 

master$window_start <-  as.numeric(master$window_start)
master$frac_LTR <-  as.numeric(master$frac_LTR)

p <- ggplot(master, aes(x = window_start, y = frac_LTR, group = chrom)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~ chrom, ncol = 1, scales = "free") +
  labs(x = "Window start", y = "Fraction of LTR") +
  theme_bw()

ggsave("frac_LTR_by_chromosome_2Mbp.pdf", plot = p, width = 10, height = 18)
ggsave("frac_LTR_by_chromosome_2Mbp.png", plot = p, width = 10, height = 18, dpi = 300)

#-----------------------------
p <- ggplot(
  master,
  aes(x = window_start, y = frac_LTR, color = chrom, group = chrom)
) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  labs(
    x = "Window start",
    y = "Fraction of LTR",
    color = "Chromosome"
  ) +
  theme_bw()

ggsave("frac_LTR_by_allchromosome_2Mbp.pdf", plot = p, width = 10, height = 6)
ggsave("frac_LTR_by_allchromosome_2Mbp.png", plot = p, width = 10, height = 6, dpi = 300)

#################################
# Plot TIR along genome

master$window_start <-  as.numeric(master$window_start)
master$frac_TIR <-  as.numeric(master$frac_TIR)

p <- ggplot(master, aes(x = window_start, y = frac_TIR, group = chrom)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~ chrom, ncol = 1, scales = "free") +
  labs(x = "Window start", y = "Fraction of LTR") +
  theme_bw()

ggsave("frac_TIR_by_chromosome_2Mbp.pdf", plot = p, width = 10, height = 18)
ggsave("frac_TIR_by_chromosome_2Mbp.png", plot = p, width = 10, height = 18, dpi = 300)

#-----------------------------
p <- ggplot(
  master,
  aes(x = window_start, y = frac_TIR, color = chrom, group = chrom)
) +
  geom_line(alpha = 0.7, na.rm = TRUE) +
  labs(
    x = "Window start",
    y = "Fraction TIR",
    color = "Chromosome"
  ) +
  theme_bw()

ggsave("frac_TIR_by_allchromosome_2Mbp.pdf", plot = p, width = 10, height = 6)
ggsave("frac_TIR_by_allchromosome_2Mbp.png", plot = p, width = 10, height = 6, dpi = 300)

#################################
# Plot SCE count along genome

master$window_start <-  as.numeric(master$window_start)
master$sce_count <-  as.numeric(master$sce_count)

p <- ggplot(master, aes(x = window_start, y = sce_count, group = chrom)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~ chrom, ncol = 1, scales = "free") +
  labs(x = "Window start", y = "SCE count") +
  theme_bw()

ggsave("SCE_by_chromosome_2Mbp.pdf", plot = p, width = 10, height = 18)
ggsave("SCE_by_chromosome_2Mbp.png", plot = p, width = 10, height = 18, dpi = 300)


#################################
# Plot AB compartments along genome

master$window_start <-  as.numeric(master$window_start)
master$comp_val <-  as.numeric(master$comp_val)

p <- ggplot(master, aes(x = window_start, y = comp_val, group = chrom)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~ chrom, ncol = 1, scales = "free") +
  labs(x = "Window start", y = "A/B compartment value") +
  theme_bw()

ggsave("AB0.5_by_chromosome_2Mbp.pdf", plot = p, width = 10, height = 18)
ggsave("AB0.5_by_chromosome_2Mbp.png", plot = p, width = 10, height = 18, dpi = 300)

######################################
# Plot correlation A/B comp value with each repeat type

# 1) Faceted boxplots: each factor split by SCE present/absent
factor_cols <- grep("^frac_", names(master), value = TRUE)

factor_cols
 # [1] "frac_CACTA_TIR_transposon"         "frac_CRE_LINE_retrotransposon"
 # [3] "frac_Copia_LTR_retrotransposon"    "frac_DNA_transposon"
 # [5] "frac_Gypsy_LTR_retrotransposon"    "frac_L1_LINE_retrotransposon"
 # [7] "frac_LTR_retrotransposon"          "frac_Mutator_TIR_transposon"
 # [9] "frac_Other"                        "frac_PIF_Harbinger_TIR_transposon"
# [11] "frac_POLE_TIR_transposon"          "frac_RTE_LINE_retrotransposon"
# [13] "frac_hAT_TIR_transposon"           "frac_helitron"
# [15] "frac_pararetrovirus"               "frac_repeat_fragment"
# [17] "frac_gene"                         "frac_mean_methylation"


# Save plots
out_base <- "repeat_ABcomp"

# --- Helper: save both png and pdf ---
save_both <- function(p, base, w=14, h=8, dpi=300) {
  ggsave(paste0(base, "_2Mbp.png"), plot = p, width = w, height = h, dpi = dpi)
  ggsave(paste0(base, "_2Mbp.pdf"), plot = p, width = w, height = h)
}


plot_data <- master %>%
  pivot_longer(
    cols = all_of(factor_cols),
    names_to = "factor",
    values_to = "value"
  ) %>%
  # group_by(factor) %>%
  # filter(
    # value >= quantile(value, 0.01, na.rm = TRUE),
    # value <= quantile(value, 0.99, na.rm = TRUE)
  # ) %>%
  ungroup()

p_scatter <- ggplot(plot_data, aes(x = comp_val, y = value)) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
  facet_wrap(~ factor, scales = "free_y", ncol = 4) +
  theme_bw() +
  labs(x = "A/B compartment", y = NULL)
  
save_both(p_scatter, paste0(out_base, "_scatter_byAB"))


################################################