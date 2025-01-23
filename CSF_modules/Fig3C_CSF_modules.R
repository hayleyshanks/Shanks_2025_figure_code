# Code for CSF proteomic module analysis, LM11A-31 Phase 2a follow-up manuscript
# Figure 3C - Placebo vs Drug comparisons of Module changes
# Written by Hayley R.C. Shanks. Last update: Jan 23, 2025
################################################################################
# SET PATHS AND LOAD PACKAGES
# install pacman if needed
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
} 
# load required libraries with pacman
pacman::p_load(rio, dplyr, ggpubr, tidyverse)
# paths to directories, code, and datasheets
wd <- "D:/HS3_ptx/code/R_codes/LM11A31_followup_paper/publication_quality/"
plt_src_code <- "csf_module_template_plot.R"
# where the final figures should be output
fig_path <- "D:/HS3_ptx/figs_graphs/Jama_neurol_2024/csf_proteomics/module_responses/"
# where the spreadsheets containing the raw data are stored
data_path <- "D:/HS3_ptx/spreadsheets/emtherapro_originals"
# spreadsheets containing the mapping of samples to subjects & the module data itself
trait_file <- "2023.12.15_PTX_Traits_Expanded_v1.5.csv"
module_file <- "PTX_CSF_ModuleEigenGenes_ALL_Drug-Placebo_Protein_Significant_Threshold0.05.csv"
################################################################################
setwd(wd)
source(plt_src_code)
meta <- import(paste(data_path, trait_file, sep="/"))
df <- import(paste(data_path, module_file, sep="/"))
###############################################################################
# FORMAT DATAFRAME

# Rename modules by module number, not by colour
df <- df %>% rename("M1"="turquoise", "M2"="blue", "M3"="brown", 
                    "M4"="yellow", "M5"="green", "M6"="red","M7"="black", 
                    "M8"="pink", "M9"="magenta", "M10"="purple")
# merge the module df with the file containing sample ID mapping
df <- left_join(df, meta[, c("STD_PRIMARYKEY", "STD_PATIENT_ID", "STD_SAMPLE_ID", "STD_TREATMENT")], by="STD_PRIMARYKEY")
# add a column for pooled drug group (placebo vs drug, no dose)
df[["Group"]] <- "drug"
df[df[["STD_TREATMENT"]] == "Placebo", "Group"] <- "placebo"
df[["Group"]] <- factor(df[["Group"]], levels = c("placebo", "drug"))
# extract visit info from the sample ID
df[["visit"]] <- sub('.*_', '', df[["STD_SAMPLE_ID"]])
# compute a difference score for V5-V1 across subjects
v1 <- df[df[["visit"]] == "V1",]
v5 <- df[df[["visit"]] == "V5",]
# we only want subjects with data for both times
to_keep <- intersect(v1[["STD_PATIENT_ID"]], v5[["STD_PATIENT_ID"]])
v1 <- v1[v1[["STD_PATIENT_ID"]] %in% to_keep,]
v5 <- v5[v5[["STD_PATIENT_ID"]] %in% to_keep,]
# sort ascending
v1 <- v1[order(v1[["STD_PATIENT_ID"]]),]
v5 <- v5[order(v5[["STD_PATIENT_ID"]]),]
###############################################################################
# CALCULATE DIFFERENCE SCORES AND PLOT RESULTS
if (all(v1[["STD_PATIENT_ID"]] == v5[["STD_PATIENT_ID"]])){
  modules <- colnames(v1)[grep("^M", colnames(df))]
  # set arguments for plot function
  to_compare <- list(c("placebo", "drug")) 
  c <- 2 # cex for ggbeeswarm
  w <- 1.5 # plot width in units
  h <- 2 # plot height in units
  font_sz <- 7 # font size for final plots
  for (m in modules){
    var_name <- paste0(m, "_difference")
    v1[[var_name]] <- v5[[m]] - v1[[m]] # difference calculation
    fig_name <- paste(fig_path, m, "_diff_wilcox_test.svg", sep = "")
    csf_module_template_plot(v1, var_name, to_compare, c, w, h, fig_name, "26-week Difference", m, font_sz)
  } 
} else{
  print("mismatched data at baseline and final visits")
}
