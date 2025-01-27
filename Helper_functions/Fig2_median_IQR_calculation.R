# Function to calculate the median and interquartile range for each of the plasma 
# variables from the LM11A-31 follow-up paper. summary data will be saved as 
# a CSV file in the data path directory
# code by Hayley R.C. Shanks, last updated Jan 27, 2024
################################################################################
# working directory
wd <- "D:/HS3_ptx/code/R_codes/LM11A31_followup_paper/publication_quality"
# path with data files from the clinical trial
data_path <- "D:/HS3_ptx/spreadsheets"
# spreadsheet with the study data
sheet_name <- "all_PTX_with_drug.csv"
# names of the plasma variables in the data sheet
vars_to_run <- c("C2N_p181_14N_APC", "C2N_p217_14N_APC", "Blennow_pTau231_APC", "Blennow_GFAP_APC")
# variable encoding participant groups
group_var <- "drugGroup"
# should outliers be excluded from calculations
outliers <- TRUE
# name of the spreadsheet to save
output_name <- "PTX_follow_up_Fig2_median_IQR_by_group.csv"
################################################################################
setwd(wd)
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(rio, dplyr)
# import the main dataset
df <- import(paste(data_path, sheet_name, sep="/"))
# encode the group_var column as a factor
df[[group_var]] <- as.factor(df[[group_var]])
# levels of the factor are our groups of interest
groups <- levels(df[[group_var]])
# empty data frame to store the output (median and IQRs)
out_df <- data.frame(Group = character(),
                     Variable = character(),
                     Median = numeric(),
                     IQR = numeric(),
                     stringsAsFactors = FALSE)
for (yvar in vars_to_run){
  if (!yvar %in% names(df)) {
    warning(paste("Variable", yvar, "not found in dataset. Skipping."))
    next
  }
  if (outliers == TRUE) {
    # nan any outliers. Outlier variables were defined using Hampel filter 
    outlier_var <- paste0(yvar, "_outliers")
    if (outlier_var %in% names(df)) {
      df[df[[outlier_var]] != 0, yvar] <- NA
    } else {
      warning(paste("Outlier variable", outlier_var, "not found. Skipping outlier handling for", yvar, "."))
    }
  }
  # calculate the median and IQR of the variable in each group
   for (group in groups) {
     vec <- df[df[[group_var]] == group, yvar]
     med <- median(vec, na.rm = TRUE)
     iq_rng <- IQR(vec, na.rm = TRUE)
     out_df <- rbind(out_df, data.frame(Group = group, Variable = yvar, Median = med, IQR = iq_rng))
  }
}
# save the out_df to file
export(out_df, file.path(data_path, output_name))
