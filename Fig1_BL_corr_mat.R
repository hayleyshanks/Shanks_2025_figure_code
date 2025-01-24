# Code for FDR-corrected correlation plot, LM11A-31 Phase 2a follow-up manuscript
# Figure 1 Baseline Correlations
# Written by Hayley R.C. Shanks. Last update: Jan 24, 2025
################################################################################
# SET PATHS AND LOAD PACKAGES
# install pacman if needed
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
} 
# load required libraries with pacman
pacman::p_load(rio, corrplot, ggplot2)
# set paths and input files
wd <- "D:/HS3_ptx/spreadsheets/"
fig_path <- "D:/HS3_ptx/figs_graphs/follow_up_publication/corr_plot/"
# original correlation matrix was created and written out in MATLAB from Shanks et al.
# 2024 Nat Med. Here we are extending and formatting the original
rho_file <- "PTX_baseline_rho_with_plasma_and_comp_and_emtherapro_subset_for_paper2.csv"
pval_file <- "PTX_baseline_p_with_plasma_and_comp_and_emtherapro_subset_for_paper2.csv"
output_name <- "baseline_correlations_lower_half_matrix.pdf"
# will be adding rectangles onto correlation matrix to group variables by data type (e.g. 
# plasma, cognitive.) set the values for the edges of the squares below
edge_vars <- c("MMSE","MEM", "MRI", "p-tau181", "M1", "AB40", "AChE")
################################################################################
setwd(wd)
rho <- import(rho_file)
p <- import(pval_file)
rownames(rho) <- colnames(rho)
p <- as.matrix(p)
# FDR adjust for multiple comparisons
p_adjusted <- p.adjust(p, method = "fdr")
# reshape vector back to original size
p_adjusted <- matrix(p_adjusted,nrow = nrow(p),ncol = ncol(p))
p_adjusted <- as.matrix(p_adjusted)
# set the p values of the diagonal to 0 so it shows up red.
diag(p_adjusted) <- 0
# carry over variable names from original spreadsheets
rownames(p_adjusted) <- rownames(rho)
colnames(p_adjusted) <- colnames(rho)
# blue to red gradient colour palette
cols <- colorRampPalette(c("blue", "white", "red"))(200)
# save output as a pdf
pdf(paste0(fig_path, output_name), width=8, height =9)
corrplot(as.matrix(rho), type = "lower", diag = TRUE, method = "color", title = "", 
         p.mat = p_adjusted, sig.level = 0.05, insig = "blank", col = cols, tl.pos = "lt",  
         tl.col = "black", tl.cex = 0.8, addgrid.col = 'grey', cl.cex = 0.8) |>
  corrRect(name = edge_vars, lwd =1)
dev.off()
