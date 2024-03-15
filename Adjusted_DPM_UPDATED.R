install.packages("rstan", type = "source")
install.packages("ggforce")
install.packages("ggrepel")
install.packages("openxlsx")


library(rstan)
library(readr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(openxlsx)
options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)


GLP1Persistency_FM <- read_csv("GLP1Persistency_sampled_data.csv")
View(GLP1Persistency_FM)
print(mean(GLP1Persistency_FM$MemberAge))
print(mean(GLP1Persistency_FM$DaysToFirstMiss))


##### Selecting columns for the model
data_matrix <- as.matrix(GLP1Persistency_FM[, c("DaysToFirstMiss", "MemberAge")])  # Adjust column names as needed

#### Prepare your data for Stan
stan_data <- list(N = nrow(data_matrix), D = ncol(data_matrix), K = 20, y = data_matrix)  # Adjust K as necessary

##### Run the Stan model
fit <- stan(
  file = "2Adjusted_DPM.stan",  # Specify the path to your Stan model file
  data = stan_data,
  chains = 4,
  iter = 4000  # Increased max_treedepth to avoid warning
)

##### Check diagnostics
print(fit)
check_hmc_diagnostics(fit)
plot(fit)

####Extract fits to data frames####
posterior_samples <- extract(fit)
cluster_assignments <- posterior_samples$cluster_assignments


write.xlsx(final_cluster_assignments, file = "Final_20_cluster_assignments.xlsx", rowNames=TRUE)


####CLuster Assignment Probability Analysis######
install.packages("pheatmap")
library(pheatmap)
library(tidyr)
library(ggplot2)
cluster20_prob<-read_csv("Cluste_20_only_proba.csv")
cluster13_prob<-read_csv("Cluster13_only_Prob.csv")


# Extracting patient IDs
patient_ids <- cluster20_prob$`Patient ID`

# Creating a matrix for the heatmap, excluding the patient ID column
prob_matrix <- as.matrix(cluster20_prob[,-1])
prob_matrix2 <- as.matrix(cluster20_prob_sub[,-1])

# Setting row names of the matrix to patient IDs
rownames(prob_matrix) <- patient_ids
rownames(prob_matrix2)< patient_ids2

#####Heatmap Probabilities#########
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
install.packages("clue")
library(clue)
library(circlize)

# Assuming prob_matrix is your probability matrix
# Transpose it if necessary to align individuals on the y-axis and clusters on the x-axis
prob_matrix <- t(prob_matrix)
color_gradient <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

# Create the heatmap
Heatmap(prob_matrix,
        name = "Probability",
        col = color_gradient,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_rot = 45, # Rotate cluster names for better visibility
        row_names_side = "left", # Patient IDs on the left
        column_title = "Cluster Probability",
        row_title = "Patient ID")
# Generate the heatmap Subsets

prob_matrix2_subset <- prob_matrix[1:30, ]

Heatmap(prob_matrix2_subset,
        name = "Probability",
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), # Customize colors as needed
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_names_side = "top", # Display cluster IDs at the top
        row_names_side = "left", # Display individual IDs on the left
        row_title = "Individuals",
        column_title = "Clusters",
        heatmap_legend_param = list(title = "Membership Probability", at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"))
)



