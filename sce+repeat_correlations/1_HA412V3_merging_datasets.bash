###################################
# Merging datasets in R for repeat, SCE and A/B compartment correlations
# November 2026
####################################

# Averages for 2 Mbp windows: A/B compartments (0.5Mbp -> 2Mbp), 
# methylation (100kbp->2Mbp), and  short read counts (1Mbp -> 2Mbp)
# Calculates per-class repeat and gene coverage fractions for each 2Mbp window
# Counts and presence/absence of SCE breakpoints in each 2 Mbp window
# Masks telomeric regions based on mask bed file - values within bed file range are set to NA

#############################################

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

# ----------------------------
# File paths
# ----------------------------
f_ab            <- "Ha412_ABcompval_0.5Mbp.txt"
f_repeat_annot  <- "TEsorter_annotation.ALLchrs_categorized.txt"
f_genes         <- "Ha412H0V3_genes.bed"
f_meth          <- "HA412HO.5mC.windows_100000.circos.bed"
f_counts        <- "strandseq_readcounts_in_1mb_bins_v3.bed"
f_sce           <- "breakPointSummary_switchesMin5Mb_verygoodqualonly.txt"
f_mask          <- "telomere_mask_for_5Mb_bins.bed"

win_size <- 2000000
fout <- "master_final.csv"

# ----------------------------
# functions
# ----------------------------
norm_chrom <- function(ch) {
  ch <- as.character(ch)
  m <- sub(".*(Chr0*([0-9]+)).*$", "Chr\\2", ch)
  ifelse(grepl("^Chr[0-9]+$", m), m, NA_character_)
}

read_tab <- function(path, header = FALSE) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!grepl("^\\s*#", lines)]  # discard comment lines
  fread(
    text = paste(lines, collapse = "\n"),
    header = header,
    sep = "\t",
    fill = TRUE
  )
}

# ----------------------------
# Read mask file and merge mask intervals
# ----------------------------
mask <- read_tab(f_mask, header = FALSE)
if (ncol(mask) < 3) stop("Mask file has fewer than 3 columns.")

setnames(mask, 1:3, c("chrom_raw", "start0", "end0"))
mask[, `:=`(
  chrom = norm_chrom(chrom_raw),
  start0 = as.numeric(start0),
  end0 = as.numeric(end0)
)]
mask <- mask[!is.na(chrom) & !is.na(start0) & !is.na(end0) & end0 > start0]

# Derive lengths from the maximum mask end per chromosome
chr_lengths <- mask[
  , .(chr_len = max(end0, na.rm = TRUE)),
  by = chrom
]

# Clip masks to chromosome bounds
mask <- merge(mask, chr_lengths, by = "chrom", all = FALSE)
mask[, `:=`(start0 = pmax(0, start0), end0 = pmin(chr_len, end0))]
mask <- mask[end0 > start0]

# Merge overlapping or adjacent mask intervals per chromosome
mask <- mask[order(chrom, start0, end0)]
mask[, prior_max_end := shift(cummax(end0), fill = -Inf), by = chrom]
mask[, mask_group := cumsum(start0 > prior_max_end), by = chrom]
mask_merged <- mask[
  , .(start0 = min(start0), end0 = max(end0)),
  by = .(chrom, mask_group)
][, mask_group := NULL]

# ----------------------------
# Build 2 Mbp windows and determine unmasked bases
# ----------------------------

backbone <- chr_lengths[
  , .(window_index = seq_len(ceiling(chr_len / win_size))),
  by = .(chrom, chr_len)
]

backbone[, window_start := (window_index - 1) * win_size]
backbone[, window_end := pmin(window_index * win_size, chr_len)]
backbone[, window_bp := window_end - window_start]
backbone[, genomic_mid := (window_start + window_end) / 2]

# Intersect merged masks with windows.
mask_exp <- merge(
  mask_merged,
  backbone[, .(chrom, window_index, window_start, window_end)],
  by = "chrom",
  allow.cartesian = TRUE
)
mask_exp[, masked_bp := pmax(
  0,
  pmin(end0, window_end) - pmax(start0, window_start)
)]
masked_win <- mask_exp[
  , .(masked_bp = sum(masked_bp)),
  by = .(chrom, window_index)
]

backbone <- merge(
  backbone, masked_win,
  by = c("chrom", "window_index"),
  all.x = TRUE
)
backbone[is.na(masked_bp), masked_bp := 0]
backbone[, unmasked_bp := window_bp - masked_bp]

# ----------------------------
# Read A/B compartments 
# aggregate to 2 Mbp windows excluding masked bases
# ----------------------------

ab <- read_tab(f_ab, header = FALSE)

if (ncol(ab) < 5) {
  stop("Compartment file has fewer than 5 columns.")
}

setnames(ab, 1:5, c(
  "chrom_raw", "start", "end", "compartment", "comp_val"
))

ab[, `:=`(
  chrom = norm_chrom(chrom_raw),
  start = as.numeric(start),
  end = as.numeric(end),
  compartment = as.character(compartment),
  comp_val = as.numeric(comp_val)
)]

ab <- ab[
  !is.na(chrom) &
    !is.na(start) &
    !is.na(end) &
    !is.na(comp_val) &
    end >= start
]

# Convert 1-based inclusive compartment coordinates to
# 0-based, half-open coordinates.
ab[, `:=`(
  start0 = start - 1,
  end0 = end
)]

# Find the 2 Mb windows each interval overlaps.
ab[, `:=`(
  win_start = floor(start0 / win_size) + 1,
  win_end = floor((end0 - 1) / win_size) + 1
)]

# Expand intervals across overlapping windows.
ab_exp <- ab[
  ,
  .(window_index = seq.int(win_start[1], win_end[1])),
  by = .(
    chrom, start0, end0, compartment, comp_val,
    win_start, win_end
  )
]

# Get the actual window boundaries and calculate overlap bases.
ab_exp <- merge(
  ab_exp,
  backbone[, .(chrom, window_index, window_start, window_end)],
  by = c("chrom", "window_index"),
  all = FALSE
)

ab_exp[, `:=`(
  seg_start = pmax(start0, window_start),
  seg_end = pmin(end0, window_end)
)]
ab_exp[, overlap_bp := pmax(0, seg_end - seg_start)]
ab_exp[, ab_row := .I]

# Find the overlap between each compartment segment and the mask.
ab_mask <- merge(
  ab_exp[, .(ab_row, chrom, seg_start, seg_end)],
  mask_merged[, .(chrom, mask_start = start0, mask_end = end0)],
  by = "chrom",
  allow.cartesian = TRUE
)

ab_mask[, masked_overlap := pmax(
  0,
  pmin(seg_end, mask_end) - pmax(seg_start, mask_start)
)]

masked_ab <- ab_mask[
  , .(masked_bp = sum(masked_overlap)),
  by = ab_row
]

# Subtract masked bases from the segment's weighting bases.
ab_exp <- merge(ab_exp, masked_ab, by = "ab_row", all.x = TRUE)
ab_exp[is.na(masked_bp), masked_bp := 0]
ab_exp[, unmasked_overlap_bp := pmax(0, overlap_bp - masked_bp)]

# Numeric compartment value: overlap-weighted mean over unmasked bases.
ab_win <- ab_exp[
  unmasked_overlap_bp > 0,
  .(
    comp_val = weighted.mean(
      comp_val,
      unmasked_overlap_bp,
      na.rm = TRUE
    )
  ),
  by = .(chrom, window_index)
]

# Optional A/B label: category covering the most unmasked bases.
ab_group_bp <- ab_exp[
  unmasked_overlap_bp > 0,
  .(unmasked_bp = sum(unmasked_overlap_bp)),
  by = .(chrom, window_index, compartment)
]

ab_label <- ab_group_bp[
  order(-unmasked_bp),
  .SD[1],
  by = .(chrom, window_index)
][
  , .(chrom, window_index, compartment)
]

ab_win <- merge(
  ab_win,
  ab_label,
  by = c("chrom", "window_index"),
  all = TRUE
)

# ----------------------------
# Repeat annotations
# Input coordinates are 1-based
# Masked bases are excluded from repeat bp and denominator
# ----------------------------
# repann <- fread(f_repeat_annot, header = TRUE, sep = "\t")
# if (ncol(repann) < 5) stop("Repeat annotation file has fewer than 5 columns.")
# setnames(repann, 1:10, c(
  # "TEName", "chrom_raw", "start", "end", "rep_class",
  # "Order", "Superfamily", "Clade", "Complete", "combined_class"
# ))

# for combined repeat annotations
repann <- fread(f_repeat_annot, header = TRUE, sep = "\t")
if (ncol(repann) < 5) stop("Repeat annotation file has fewer than 5 columns.")
setnames(repann, 1:10, c(
  "TEName", "chrom_raw", "start", "end", "EDTA_class",
  "Order", "Superfamily", "Clade", "Complete", "rep_class"
))

repann[, `:=`(
  chrom = norm_chrom(chrom_raw),
  start0 = as.numeric(start) - 1,
  end0 = as.numeric(end)
)]
repann <- repann[
  !is.na(chrom) & !is.na(start0) & !is.na(end0) & end0 > start0
]
repann[is.na(rep_class) | rep_class == "", rep_class := "Unknown"]

repann[, `:=`(
  win_start = floor(start0 / win_size) + 1,
  win_end = floor((end0 - 1) / win_size) + 1
)]
repann_exp <- repann[
  , .(window_index = seq.int(win_start, win_end)),
  by = .(chrom, rep_class, start0, end0)
]
repann_exp <- merge(
  repann_exp,
  backbone[, .(chrom, window_index, window_start, window_end)],
  by = c("chrom", "window_index"),
  all = FALSE
)
repann_exp[, `:=`(
  int_start = pmax(start0, window_start),
  int_end = pmin(end0, window_end)
)]

# Find the masked overlap for each repeat segment within its window.
te_mask <- merge(
  repann_exp[, .(chrom, window_index, rep_class, int_start, int_end)],
  mask_merged[, .(chrom, mask_start = start0, mask_end = end0)],
  by = "chrom",
  allow.cartesian = TRUE
)
te_mask[, masked_overlap := pmax(
  0, pmin(int_end, mask_end) - pmax(int_start, mask_start)
)]
masked_te <- te_mask[
  , .(masked_overlap = sum(masked_overlap)),
  by = .(chrom, window_index, rep_class, int_start, int_end)
]

repann_exp <- merge(
  repann_exp, masked_te,
  by = c("chrom", "window_index", "rep_class", "int_start", "int_end"),
  all.x = TRUE
)
repann_exp[is.na(masked_overlap), masked_overlap := 0]
repann_exp[, repeat_bp := pmax(
  0, int_end - int_start - masked_overlap
)]

# Total repeat bp and repeat fraction.
rep_cov_win <- repann_exp[
  , .(repeat_bp = sum(repeat_bp, na.rm = TRUE)),
  by = .(chrom, window_index)
]
rep_cov_win <- merge(
  rep_cov_win,
  backbone[, .(chrom, window_index, unmasked_bp)],
  by = c("chrom", "window_index"),
  all.x = TRUE
)
rep_cov_win[, repeat_frac := fifelse(
  unmasked_bp > 0, repeat_bp / unmasked_bp, NA_real_
)]

# Per-class repeat fractions. Keep the 15 most abundant classes; pool others.
rep_class_win <- repann_exp[
  , .(repeat_bp = sum(repeat_bp, na.rm = TRUE)),
  by = .(chrom, window_index, rep_class)
]
top_classes <- rep_class_win[
  , .(genome_repeat_bp = sum(repeat_bp, na.rm = TRUE)),
  by = rep_class
][order(-genome_repeat_bp)][1:min(15, .N), rep_class]

rep_class_win[, rep_class2 := fifelse(
  rep_class %in% top_classes, rep_class, "Other"
)]
rep_class_agg <- rep_class_win[
  , .(repeat_bp = sum(repeat_bp, na.rm = TRUE)),
  by = .(chrom, window_index, rep_class2)
]
rep_class_wide <- dcast(
  rep_class_agg,
  chrom + window_index ~ rep_class2,
  value.var = "repeat_bp",
  fill = 0
)
class_bp_cols <- setdiff(names(rep_class_wide), c("chrom", "window_index"))
rep_class_wide <- merge(
  rep_class_wide,
  backbone[, .(chrom, window_index, unmasked_bp)],
  by = c("chrom", "window_index"),
  all.x = TRUE
)
for (cc in class_bp_cols) {
  set(rep_class_wide, j = cc, value = fifelse(
    rep_class_wide$unmasked_bp > 0,
    rep_class_wide[[cc]] / rep_class_wide$unmasked_bp,
    NA_real_
  ))
}
setnames(rep_class_wide, class_bp_cols, paste0("frac_", class_bp_cols))
rep_class_wide[, unmasked_bp := NULL]

# ----------------------------
# SCE breakpoints
# Each SCE is assigned once, to the window containing its midpoint.
# ----------------------------
sce <- fread(f_sce, header = FALSE)

if (ncol(sce) < 7) stop("SCE BED file has fewer than 7 columns.")
setnames(sce, 1:7, c("seqnames", "start", "end", "CI.start", "CI.end", "genoT", "filenames"))


if (!all(c("seqnames", "start", "end") %in% names(sce))) {
  stop("SCE file must have seqnames, start, and end columns.")
}

sce[, `:=`(
  chrom = norm_chrom(seqnames),
  start = as.numeric(start),
  end = as.numeric(end)
)]
sce <- sce[!is.na(chrom) & !is.na(start) & !is.na(end)]
sce[, midpoint0 := ((start + end) / 2) - 1]

sce[, sce_id := .I]

sce_in_mask <- merge(
  sce[, .(sce_id, chrom, midpoint0)],
  mask_merged,
  by = "chrom",
  allow.cartesian = TRUE
)[midpoint0 >= start0 & midpoint0 < end0, unique(sce_id)]

sce <- sce[!sce_id %in% sce_in_mask]
sce[, window_index := floor(midpoint0 / win_size) + 1]

sce_counts <- sce[
  , .(sce_count = .N),
  by = .(chrom, window_index)
]
sce_counts[, sce_present := as.integer(sce_count > 0)]

# ----------------------------
# Genes
# Input is a .bed file: 0-based, half-open
# ----------------------------
genes <- read_tab(f_genes, header = FALSE)
if (ncol(genes) < 4) stop("Gene BED file has fewer than 4 columns.")
setnames(genes, 1:4, c("chrom_raw", "start0", "end0", "geneID"))
genes[, `:=`(
  chrom = norm_chrom(chrom_raw),
  start0 = as.numeric(start0),
  end0 = as.numeric(end0)
)]
genes <- genes[!is.na(chrom) & end0 > start0]
genes[, `:=`(
  win_start = floor(start0 / win_size) + 1,
  win_end = floor((end0 - 1) / win_size) + 1
)]
genes_exp <- genes[
  , .(window_index = seq.int(win_start, win_end)),
  by = .(chrom, geneID, start0, end0)
]
genes_exp <- merge(
  genes_exp,
  backbone[, .(chrom, window_index, window_start, window_end)],
  by = c("chrom", "window_index"),
  all = FALSE
)
genes_exp[, `:=`(
  int_start = pmax(start0, window_start),
  int_end = pmin(end0, window_end)
)]
gene_win <- genes_exp[
  , .(
    gene_count = uniqueN(geneID),
    gene_bp = sum(pmax(0, int_end - int_start), na.rm = TRUE)
  ),
  by = .(chrom, window_index)
]
gene_win <- merge(
  gene_win,
  backbone[, .(chrom, window_index, unmasked_bp)],
  by = c("chrom", "window_index"),
  all.x = TRUE
)
gene_win[, frac_gene := fifelse(
  unmasked_bp > 0, gene_bp / unmasked_bp, NA_real_
)]
gene_win[, unmasked_bp := NULL]

# ----------------------------
# Methylation: average 100 kb window values by 2 Mbp window
# Input coordinates are in a .bed file, 0-based, half-open
# ----------------------------
meth <- read_tab(f_meth, header = FALSE)
if (ncol(meth) < 4) stop("Methylation file has fewer than 4 columns.")
setnames(meth, 1:4, c("chrom_raw", "start0", "end0", "methylation_pct"))
meth[, `:=`(
  chrom = norm_chrom(chrom_raw),
  start0 = as.numeric(start0),
  end0 = as.numeric(end0),
  methylation_pct = as.numeric(methylation_pct)
)]
meth <- meth[!is.na(chrom) & !is.na(methylation_pct) & end0 > start0]
meth[, `:=`(
  win_start = floor(start0 / win_size) + 1,
  win_end = floor((end0 - 1) / win_size) + 1
)]
meth_exp <- meth[
  , .(window_index = seq.int(win_start, win_end)),
  by = .(chrom, start0, end0, methylation_pct)
]
meth_exp <- merge(
  meth_exp,
  backbone[, .(chrom, window_index, window_start, window_end)],
  by = c("chrom", "window_index"),
  all = FALSE
)
meth_exp[, overlap_bp := pmax(
  0, pmin(end0, window_end) - pmax(start0, window_start)
)]
meth_win <- meth_exp[
  , .(
    frac_mean_methylation = weighted.mean(
      methylation_pct, overlap_bp, na.rm = TRUE
    ),
    methylation_n_100kb = .N
  ),
  by = .(chrom, window_index)
]

# ----------------------------
# Short-read counts
# The bins are assigned by start coordinate
# ----------------------------
counts <- read_tab(f_counts, header = FALSE)
if (ncol(counts) < 4) stop("Read-count file has fewer than 4 columns.")
setnames(counts, 1:4, c("chrom_raw", "start0", "end0", "read_count"))
counts[, `:=`(
  chrom = norm_chrom(chrom_raw),
  start0 = as.numeric(start0),
  end0 = as.numeric(end0),
  read_count = as.numeric(read_count)
)]
counts <- counts[!is.na(chrom) & !is.na(read_count)]
counts[, window_index := floor(start0 / win_size) + 1]
counts_win <- counts[
  , .(read_count_sum = sum(read_count, na.rm = TRUE)),
  by = .(chrom, window_index)
]

# ----------------------------
# Join datasets
# ----------------------------
master <- copy(backbone)

master <- merge(
  master, ab_win,
  by = c("chrom", "window_index"), all.x = TRUE
)
master <- merge(
  master, rep_cov_win[, .(chrom, window_index, repeat_bp, repeat_frac)],
  by = c("chrom", "window_index"), all.x = TRUE
)
master <- merge(
  master, rep_class_wide,
  by = c("chrom", "window_index"), all.x = TRUE
)
master <- merge(
  master, sce_counts,
  by = c("chrom", "window_index"), all.x = TRUE
)
master <- merge(
  master, gene_win,
  by = c("chrom", "window_index"), all.x = TRUE
)
master <- merge(
  master, meth_win,
  by = c("chrom", "window_index"), all.x = TRUE
)
master <- merge(
  master, counts_win,
  by = c("chrom", "window_index"), all.x = TRUE
)

master[is.na(sce_count), sce_count := 0L]
master[is.na(sce_present), sce_present := 0L]

master[unmasked_bp == 0, `:=`(
  sce_count = NA_integer_,
  sce_present = NA_integer_
)]

# A fully masked window has no usable denominator.
master[unmasked_bp <= 0, `:=`(
  repeat_frac = NA_real_,
  repeat_frac_total = NA_real_
)]
master[, repeat_frac_total := repeat_frac]


# After all datasets have been joined into master - label masked regions:
master[, mask_touched := !is.na(masked_bp) & masked_bp > 0]
masked_rows <- which(master$mask_touched)

dataset_cols <- setdiff(
  names(master),
  c(
    "chrom", "window_index", "chr_len",
    "window_start", "window_end", "window_bp", "genomic_mid",
    "masked_bp", "unmasked_bp", "mask_touched"
  )
)

for (cc in dataset_cols) {
  set(master, i = masked_rows, j = cc, value = NA)
}


# Save the complete joined table
master <- master[!is.na(chrom) & chrom != "Chr0"]
fwrite(master, fout)

# checks
print(master[, .N, by = chrom][order(chrom)])

print(master[, .(
  windows = .N,
  fully_masked_windows = sum(unmasked_bp <= 0),
  min_unmasked_bp = min(unmasked_bp, na.rm = TRUE),
  max_unmasked_bp = max(unmasked_bp, na.rm = TRUE)
)])

# check masked windows are NA not 0 SCE
master[unmasked_bp == 0,
       .(chrom, window_index, unmasked_bp, repeat_frac, repeat_frac_total)]

# check column names
names(master)
 # [1] "chrom"                             "window_index"
 # [3] "chr_len"                           "window_start"
 # [5] "window_end"                        "window_bp"
 # [7] "genomic_mid"                       "masked_bp"
 # [9] "unmasked_bp"                       "comp_val"
# [11] "compartment"                       "repeat_bp"
# [13] "repeat_frac"                       "frac_CACTA_TIR_transposon"
# [15] "frac_CRE_LINE_retrotransposon"     "frac_Copia_LTR_retrotransposon"
# [17] "frac_DNA_transposon"               "frac_Gypsy_LTR_retrotransposon"
# [19] "frac_L1_LINE_retrotransposon"      "frac_LTR_retrotransposon"
# [21] "frac_Mutator_TIR_transposon"       "frac_Other"
# [23] "frac_PIF_Harbinger_TIR_transposon" "frac_POLE_TIR_transposon"
# [25] "frac_RTE_LINE_retrotransposon"     "frac_hAT_TIR_transposon"
# [27] "frac_helitron"                     "frac_pararetrovirus"
# [29] "frac_repeat_fragment"              "sce_count"
# [31] "sce_present"                       "gene_count"
# [33] "gene_bp"                           "frac_gene"
# [35] "frac_mean_methylation"             "methylation_n_100kb"
# [37] "read_count_sum"                    "repeat_frac_total"


######################################
