source("run_qstfst.r")

library(gaston)
library(hierfstat)
library(JGTeach)

args <- commandArgs(trailingOnly = TRUE)

# Set parameters
replicate_number <- as.integer(args[1])
Population_structure <- args[2]
number_of_pop <- as.integer(args[3])
generations <- as.integer(args[4])
correlation <- ifelse(args[5] == "NA", "", paste0("_", args[5]))
wdiff <- args[6]
wvar <- args[7]
experiment_name <- args[8]

# Load prepared data
folder_name <- paste0("wdiff", wdiff, "_wvar", wvar)
intermediate_dir <- paste0("intermediate/", experiment_name, "/", folder_name, "/")
intermediate_file <- paste0(intermediate_dir, "prepared_data_rep", replicate_number, ".RData")

cat("Loading data from:", intermediate_file, "\n")
load(intermediate_file)

# Run QST-FST analysis
cat("Running QST-FST analysis...\n")

qstfst_result <- run_qstfst(dataf = dataf, 
                           fst_neutral = fst.founders, 
                           ns = n_sire, 
                           nd = n_dam, 
                           no = no, 
                           np = np)

qstfst_result_Neutral <- run_qstfst(dataf = dataf_neutral, 
                                   fst_neutral = fst.founders, 
                                   ns = n_sire, 
                                   nd = n_dam, 
                                   no = no, 
                                   np = np)

# Save QST-FST results
qstfst_results <- data.frame(
  replicate_number = replicate_number,
  folder_name = folder_name,
  p_value_QSTFST = qstfst_result$p_value,
  p_value_QSTFST_Neutral = qstfst_result_Neutral$p_value
)

# Save results
results_dir <- paste0("results/",experiment_name, "/")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

qstfst_file <- paste0(results_dir, "qstfst_results_", folder_name, ".csv")

write.table(
  qstfst_results, file = qstfst_file, append = TRUE, 
  sep = ",", col.names = !file.exists(qstfst_file), row.names = FALSE
)

cat("QST-FST analysis completed. Results saved to:", qstfst_file, "\n")