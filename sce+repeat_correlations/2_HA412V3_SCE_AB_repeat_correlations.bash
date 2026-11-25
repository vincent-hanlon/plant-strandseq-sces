###################################
# Correlations of repeat, SCE and A/B compartment along each chromosome
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


####################################
# What factors have a different distribution when SCEs are present versus absent?

# baseline model for SCE presence
m_base <- glmmTMB(sce_present ~  (1|chrom),
                  family = binomial, data = master)

repeat_cols <- grep("^frac_", names(master), value = TRUE)

repeat_cols <- c(repeat_cols, "comp_val")

# fit model: add one repeat type or gene_frac
models <- lapply(repeat_cols, function(rc) {
  f <- as.formula(paste0(
    "sce_present ~ `", rc, "` + (1|chrom)"
  ))
  glmmTMB(f, family = binomial, data = master)
})

names(models) <- repeat_cols

# LRT p-values comparing to baseline (nested models)
pvals <- rbindlist(lapply(repeat_cols, function(rc) {
  m <- models[[rc]]
  a <- anova(m_base, m)  # row 2 is the larger model
  data.table(
    repeat_type = rc,
    df_diff = a$Df[2],
    AIC_base = AIC(m_base),
    AIC_model = AIC(m),
    p_value = a$`Pr(>Chisq)`[2]
  )
}))

pvals[order(p_value)]

# combine LTR annotations, 2 Mbp windows, 5Mbp cut
                 # repeat_type df_diff AIC_base AIC_model    p_value direction
                      # <char>   <num>    <num>     <num>      <num>    <char>
# 1:                 frac_gene       3 1993.166  1988.755 0.01134473  negative
# 2:                frac_OTHER       3 1993.166  1989.072 0.01356442  negative
# 3:                  frac_LTR       3 1993.166  1989.828 0.02086765  positive
# 4:           0.5Mbp_comp_val       3 1993.166  1990.569 0.03203433  positive
# 5:     frac_mean_methylation       3 1993.166  1992.623 0.11075995  positive
# 6:                  frac_TIR       3 1993.166  1993.941 0.26832547  negative
# 7: frac_long_terminal_repeat       2 1993.166  1993.166 1.00000000      <NA>


# 2Mbp windows, 5Mbp telomere cut 
                          # repeat_type df_diff AIC_base AIC_model     p_value
                               # <char>   <num>    <num>     <num>       <num>
 # 1:                     frac_helitron       3 1993.166  1985.785 0.002192278
 # 2:                         frac_gene       3 1993.166  1988.755 0.011344729
 # 3:               frac_DNA_transposon       3 1993.166  1989.005 0.013063027
 # 4:                        frac_Other       3 1993.166  1989.486 0.017164593
 # 5:                   0.5Mbp_comp_val       3 1993.166  1990.569 0.032034327
 
 # 6:    frac_Copia_LTR_retrotransposon       3 1993.166  1991.524 0.056356873
 # 7:    frac_Gypsy_LTR_retrotransposon       3 1993.166  1992.231 0.086682194
 # 8:          frac_LTR_retrotransposon       3 1993.166  1992.452 0.099449185
 
 # 9:             frac_mean_methylation       3 1993.166  1992.623 0.110759953
# 10: frac_PIF_Harbinger_TIR_transposon       3 1993.166  1992.747 0.119868563
# 11:              frac_repeat_fragment       3 1993.166  1992.809 0.124720597
# 12:           frac_hAT_TIR_transposon       3 1993.166  1992.836 0.126941724
# 13:               frac_pararetrovirus       3 1993.166  1993.213 0.162244295
# 14:     frac_RTE_LINE_retrotransposon       3 1993.166  1993.333 0.175824002
# 15:     frac_CRE_LINE_retrotransposon       3 1993.166  1993.868 0.254616402
# 16:      frac_L1_LINE_retrotransposon       3 1993.166  1994.294 0.350366814
# 17:       frac_Mutator_TIR_transposon       3 1993.166  1995.146 0.886639131
# 18:         frac_CACTA_TIR_transposon       3 1993.166        NA          NA
# 19:          frac_POLE_TIR_transposon       3 1993.166        NA          NA

######################################
# Look at specific correlations
# Is there a significant difference in the LTR frac or the AB comp val 
# for windows with SCEs present vs absent

# null model
m_null <- lme4::glmer(
  sce_present ~ 1 + (1 | chrom),
  data = master,
  family = binomial
)

#--------------
# AB compartments alone
m_comp <- lme4::glmer(
  sce_present ~ comp_val + (1 | chrom),
  data = master,
  family = binomial
)
anova(m_null, m_comp, test = "Chisq")

# 2Mbp window, 5Mbp telomere cut
# Data: master
# Models:
# m_null: sce_present ~ 1 + (1 | chrom)
# m_comp: sce_present ~ comp_val + (1 | chrom)
       # npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
# m_null    2 1993.2 2003.7 -994.58    1989.2
# m_comp    3 1990.6 2006.4 -992.28    1984.6 4.5967  1    0.03203 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#---------------------------------
# Other LTR categories

m_ltr <- lme4::glmer(
  sce_present ~ frac_LTR_retrotransposon + (1 | chrom),
  data = master,
  family = binomial
)

anova(m_null, m_ltr, test = "Chisq")

# 2Mbp window, 5Mbp telomere cut
# Data: master
# Models:
# m_null: sce_present ~ 1 + (1 | chrom)
# m_ltr: sce_present ~ frac_LTR_retrotransposon + (1 | chrom)
       # npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
# m_null    2 1993.2 2003.7 -994.58    1989.2
# m_ltr     3 1992.5 2008.3 -993.23    1986.5 2.7144  1    0.09945 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#---------------------------------
# Helitrons

m_helit <- lme4::glmer(
  sce_present ~ frac_helitron + (1 | chrom),
  data = master,
  family = binomial
)

anova(m_null, m_helit, test = "Chisq")

# 2Mbp window, 5Mbp telomere cut
# Data: master
# Models:
# m_null: sce_present ~ 1 + (1 | chrom)
# m_helit: sce_present ~ frac_helitron + (1 | chrom)
        # npar    AIC    BIC  logLik -2*log(L)  Chisq Df Pr(>Chisq)
# m_null     2 1993.2 2003.7 -994.58    1989.2
# m_helit    3 1985.8 2001.6 -989.89    1979.8 9.3812  1   0.002192 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##########################################
# Account for A/B compartment value variation (comp_val) first and then check again

# Check all repeat types 
master[, sce_present := factor(sce_present)]
master[, sce_present:= as.integer(sce_present == "1")]  

# baseline model for SCE presence
m_base <- glmmTMB(sce_present ~ comp_val + (1|chrom),
                  family = binomial, data = master)

repeat_cols <- grep("^frac_", names(master), value = TRUE)

# fit model: add one repeat type or  gene_frac
models <- lapply(repeat_cols, function(rc) {
  f <- as.formula(paste0(
    "sce_present ~ comp_val + `", rc, "` + (1|chrom)"
  ))
  glmmTMB(f, family = binomial, data = master)
})
names(models) <- repeat_cols

# LRT p-values comparing to baseline
pvals <- rbindlist(lapply(repeat_cols, function(rc) {
  m <- models[[rc]]
  a <- anova(m_base, m)  # row 2 is the larger model
  data.table(
    repeat_type = rc,
    df_diff = a$Df[2],
    AIC_base = AIC(m_base),
    AIC_model = AIC(m),
    p_value = a$`Pr(>Chisq)`[2]
  )
}))

pvals[order(p_value)]

# combined annotation, 2Mbp windows + 5Mbp mask
                 # repeat_type df_diff AIC_base AIC_model    p_value
                      # <char>   <num>    <num>     <num>      <num>
# 1:                 frac_gene       4 1990.569  1989.858 0.09961061
# 2:                frac_OTHER       4 1990.569  1989.901 0.10235843
# 3:                  frac_LTR       4 1990.569  1991.595 0.32356342
# 4:                  frac_TIR       4 1990.569  1992.542 0.86936020
# 5:     frac_mean_methylation       4 1990.569  1992.566 0.95342436
# 6: frac_long_terminal_repeat       3 1990.569  1990.569 1.00000000

# 2Mbp window, 5Mbp telomere cut
                          # repeat_type df_diff AIC_base AIC_model    p_value
                               # <char>   <num>    <num>     <num>      <num>
 # 1:                     frac_helitron       4 1990.569  1985.455 0.00764806
 # 2:                        frac_Other       4 1990.569  1987.813 0.02919117
 
 # 3:                         frac_gene       4 1990.569  1989.858 0.09961061
 # 4:    frac_Copia_LTR_retrotransposon       4 1990.569  1991.059 0.21908586
 # 5:      frac_L1_LINE_retrotransposon       4 1990.569  1991.410 0.28158338
 # 6:           frac_hAT_TIR_transposon       4 1990.569  1991.806 0.38219360
 # 7:    frac_Gypsy_LTR_retrotransposon       4 1990.569  1992.218 0.55346297
 # 8: frac_PIF_Harbinger_TIR_transposon       4 1990.569  1992.418 0.69691998
 # 9:     frac_CRE_LINE_retrotransposon       4 1990.569  1992.491 0.77892992
# 10:          frac_POLE_TIR_transposon       4 1990.569  1992.530 0.84252613
# 11:         frac_CACTA_TIR_transposon       4 1990.569  1992.539 0.86077047
# 12:          frac_LTR_retrotransposon       4 1990.569  1992.564 0.93985568
# 13:                          frac_LTR       4 1990.569  1992.564 0.93985568
# 14:              frac_repeat_fragment       4 1990.569  1992.565 0.94833131
# 15:             frac_mean_methylation       4 1990.569  1992.566 0.95342436
# 16:               frac_DNA_transposon       4 1990.569        NA         NA
# 17:       frac_Mutator_TIR_transposon       4 1990.569        NA         NA
# 18:     frac_RTE_LINE_retrotransposon       4 1990.569        NA         NA
# 19:               frac_pararetrovirus       4 1990.569        NA         NA


#-----------------------------------
# Look at direction and size of effect

get_dir <- function(rc, m) {
  tab <- summary(m)$coefficients$cond
  rn <- rownames(tab)

  # match rownames 
  hits <- rn[grepl(paste0("^", rc, "$"), rn)]

# if exact match fails
  if (length(hits) == 0) {
    hits <- rn[grepl(rc, rn, fixed = TRUE)]
  }

  if (length(hits) == 0) {
    return(list(direction = NA_character_, estimate = NA_real_))
  }

  est <- tab[hits[1], "Estimate"]
  dir <- ifelse(est > 0, "positive", "negative")
  list(direction = dir, estimate = est, matched_rows = hits)
}


pvals <- rbindlist(lapply(repeat_cols, function(rc) {
  m <- models[[rc]]
  a <- anova(m_base, m)

  co <- get_dir(rc, m)

  data.table(
    repeat_type = rc,
    df_diff = a$Df[2],
    AIC_base = AIC(m_base),
    AIC_model = AIC(m),
    p_value = a$`Pr(>Chisq)`[2],
    direction = co$direction,
    beta = co$estimate
  )
}))

pvals[order(p_value)]


# Combined repeat annotation, 2 Mbp windows, 5Mbp cut
                 # repeat_type df_diff AIC_base AIC_model    p_value direction
                      # <char>   <num>    <num>     <num>      <num>    <char>
# 1:                 frac_gene       3 1993.166  1988.755 0.01134473  negative
# 2:                frac_OTHER       3 1993.166  1989.072 0.01356442  negative
# 3:                  frac_LTR       3 1993.166  1989.828 0.02086765  positive
# 4:                  comp_val       3 1993.166  1990.569 0.03203433  positive
# 5:     frac_mean_methylation       3 1993.166  1992.623 0.11075995  positive
# 6:                  frac_TIR       3 1993.166  1993.941 0.26832547  negative
# 7: frac_long_terminal_repeat       2 1993.166  1993.166 1.00000000      <NA>
           # beta
          # <num>
# 1:  -2.88055715
# 2:  -5.64104049
# 3:   1.04771472
# 4:   2.59541370
# 5:   0.06351544
# 6: -13.83986229
# 7:           NA

# repeat annotation not combined - 2 Mbp windows, 5Mbp cut
                          # repeat_type df_diff AIC_base AIC_model     p_value
                               # <char>   <num>    <num>     <num>       <num>
 # 1:                     frac_helitron       3 1993.166  1985.785 0.002192278
 # 2:                         frac_gene       3 1993.166  1988.755 0.011344729
 # 3:               frac_DNA_transposon       3 1993.166  1989.005 0.013063027
 # 4:                        frac_Other       3 1993.166  1989.486 0.017164593
 # 5:    frac_Copia_LTR_retrotransposon       3 1993.166  1991.524 0.056356873
 # 6:    frac_Gypsy_LTR_retrotransposon       3 1993.166  1992.231 0.086682194
 # 7:          frac_LTR_retrotransposon       3 1993.166  1992.452 0.099449185
 # 8:                          frac_LTR       3 1993.166  1992.452 0.099449185
 
 
 # 9:             frac_mean_methylation       3 1993.166  1992.623 0.110759953
# 10: frac_PIF_Harbinger_TIR_transposon       3 1993.166  1992.747 0.119868563
# 11:              frac_repeat_fragment       3 1993.166  1992.809 0.124720597
# 12:           frac_hAT_TIR_transposon       3 1993.166  1992.836 0.126941724
# 13:               frac_pararetrovirus       3 1993.166  1993.213 0.162244295
# 14:     frac_RTE_LINE_retrotransposon       3 1993.166  1993.333 0.175824002
# 15:     frac_CRE_LINE_retrotransposon       3 1993.166  1993.868 0.254616402
# 16:      frac_L1_LINE_retrotransposon       3 1993.166  1994.294 0.350366814
# 17:       frac_Mutator_TIR_transposon       3 1993.166  1995.146 0.886639131
# 18:         frac_CACTA_TIR_transposon       3 1993.166        NA          NA
# 19:          frac_POLE_TIR_transposon       3 1993.166        NA          NA
    # direction          beta
       # <char>         <num>
 # 1:  negative  -59.93004849
 # 2:  negative   -2.88055715
 # 3:  negative -320.49064418
 # 4:  negative -279.17295591
 # 5:  positive    3.36847726
 # 6:  positive    1.43121904
 # 7:  positive    1.35956160
 # 8:  positive    1.35956160
 
 # 9:  positive    0.06351544
# 10:  negative  -44.31400267
# 11:  negative   -5.08977410
# 12:  negative -109.17449706
# 13:  negative  -79.38908904
# 14:  negative  -70.88750318
# 15:  negative  -51.35339948
# 16:  negative   -2.67336160
# 17:  positive    3.09378986
# 18:  negative  -26.33219569
# 19:  negative  -63.52677064


######################################
# Compare repeats and A/B compartments to SCE counts 
# with linear model (continuous SCE count)

m_comp <- lmerTest::lmer(
  sce_count ~ comp_val + (1 | chrom),
  data = master
)

stats::anova(m_comp)

# 2Mbp window, 5Mbp telomere cut
# Type III Analysis of Variance Table with Satterthwaite's method
          # Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# comp_val 0.89691 0.89691     1 950.08   1.237 0.2663


#--------------------------

m_ltr <- lmerTest::lmer(
  sce_count ~ frac_LTR + (1 | chrom),
  data = master
)

stats::anova(m_ltr)

# 2Mbp window, 5Mbp telomere cut
# Type III Analysis of Variance Table with Satterthwaite's method
                          # Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# frac_LTR_retrotransposon 0.65679 0.65679     1 754.12  0.9066 0.3413


#---------------------------
m_helit <- lmerTest::lmer(
  sce_count ~ frac_helitron + (1 | chrom),
  data = master
)

stats::anova(m_helit)

# 2Mbp window, 5Mbp mask
# Type III Analysis of Variance Table with Satterthwaite's method
              # Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)
# frac_helitron 3.0361  3.0361     1 1414.5   4.198 0.04066 *


#####################################
# Check correlation between short read mapping and repeats/SCE detection

m_read_rep <- lmerTest::lmer(
  read_count_sum ~ repeat_frac + (1 | chrom),
  data = master
)
stats::anova(m_read_rep)

# 2Mbp window
# Type III Analysis of Variance Table with Satterthwaite's method
                # Sum Sq    Mean Sq NumDF  DenDF F value    Pr(>F)
# repeat_frac 7.5107e+10 7.5107e+10     1 1427.3  3388.6 < 2.2e-16 ***


#-------------------------------
m_read_SCE <- lmerTest::lmer(
  sce_count ~ read_count_sum + (1 | chrom),
  data = master
)
stats::anova(m_read_SCE)

# 2 Mbp window
# Type III Analysis of Variance Table with Satterthwaite's method
                   # Sum Sq    Mean Sq NumDF  DenDF F value Pr(>F)
# read_count_sum 8.6067e-05 8.6067e-05     1 1273.3   1e-04 0.9913


#-------------------
m_read_rep_SCE <- lmerTest::lmer(
  sce_count ~ read_count_sum + repeat_frac + (1 | chrom),
  data = master
)
stats::anova(m_read_rep_SCE)

# 2Mbp window
# Type III Analysis of Variance Table with Satterthwaite's method
                # Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# read_count_sum 0.53101 0.53101     1  975.3  0.7323 0.3923
# repeat_frac    0.75487 0.75487     1 1135.5  1.0410 0.3078


#############################################################