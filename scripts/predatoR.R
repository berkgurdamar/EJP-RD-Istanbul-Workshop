# installation ------------------------------------------------------------

if(!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}

install_github("berkgurdamar/predatoR")

# Gene name included ------------------------------------------------------

library(predatoR)

input_df <- as.data.frame(rbind(c("3SQJ", "A", 196, "GLN", "LEU", "ALB"),
                                c("3SQJ", "A", 396, "GLU", "LYS", "ALB")))

# ca atoms approach with 7 Angstrom
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE, distance_cutoff = 7, network_approach = "ca")

# all atoms approach with 5 Angstrom
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE, distance_cutoff = 5, network_approach = "all")

# Gene name not included --------------------------------------------------

input_df <- as.data.frame(rbind(c("3SQJ", "A", 196, "GLN", "LEU"),
                                c("3SQJ", "A", 396, "GLU", "LYS")))
# ca approach with 7 Angstrom
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = FALSE, distance_cutoff = 7, network_approach = "ca")


# Partially included gene names -------------------------------------------

input_df <- as.data.frame(rbind(c("3SQJ", "A", 196, "GLN", "LEU", "ALB"),
                                c("3SQJ", "A", 396, "GLU", "LYS", "")))

# ca approach with 7 Angstrom
pred_res <- predatoR(info_df =  input_df, n_threads = 8, gene_name_info = TRUE, distance_cutoff = 7, network_approach = "ca")


# Exploratory usage -------------------------------------------------------

# Cα approach with 7.6Å Angstrom
prediction_result <- predatoR(input_df, n_threads = 8, gene_name_info = TRUE, distance_cutoff = 7.6, network_approach = "ca") 


# AlphaFold ---------------------------------------------------------------

input_df <- as.data.frame(cbind("AF-Q8NBT3-F1-model_v4", "A", 4, "LEU", "LYS", "TMEM145"))

# ca approach with 7 Angstrom
prediction_result <- predatoR(input_df, n_threads = 8, gene_name_info = TRUE, distance_cutoff = 7, network_approach = "ca", PDB_path = "D:/Downloads/") 

