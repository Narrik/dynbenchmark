# Set seed and randomly shuffle matrix
set.seed(42)
rows <- sample(nrow(countsMatrix))
countsMatrix_shuffled <- countsMatrix[rows,]
logCountsMatrix_shuffled <- logCountsMatrix[rows,]

countsMatrix_subset1 <- countsMatrix_shuffled[2089:8354,]
countsMatrix_subset2 <- countsMatrix_shuffled[c(1:2087,4072,4177:8354),]
#countsMatrix_subset2 <- countsMatrix_shuffled[c(1:2088,4177:8354),]
countsMatrix_subset3 <- countsMatrix_shuffled[c(1:4176,6265:8354),]
countsMatrix_subset4 <- countsMatrix_shuffled[1:6266,]

logCountsMatrix_subset1 <- logCountsMatrix_shuffled[2089:8354,]
logCountsMatrix_subset2 <- logCountsMatrix_shuffled[c(1:2087,4072,4177:8354),]
#logCountsMatrix_subset2 <- logCountsMatrix_shuffled[c(1:2088,4177:8354),]
logCountsMatrix_subset3 <- logCountsMatrix_shuffled[c(1:4176,6265:8354),]
logCountsMatrix_subset4 <- logCountsMatrix_shuffled[1:6266,]

dataset_subset1<- wrap_expression(counts = countsMatrix_subset1,expression = logCountsMatrix_subset1)
dataset_subset2<- wrap_expression(counts = countsMatrix_subset2,expression = logCountsMatrix_subset2)
dataset_subset3<- wrap_expression(counts = countsMatrix_subset3,expression = logCountsMatrix_subset3)
dataset_subset4<- wrap_expression(counts = countsMatrix_subset4,expression = logCountsMatrix_subset4)

methods <- c("wishbone")
for (method in methods){
  for (i in 4){
    # Model and function variable names
    ti_dr_model <- paste(method,"_dr_model_subset",i, sep="")
    ti_method <- match.fun(paste("ti_",method, sep=""))
    # Run method
    #temp_model <- infer_trajectory(eval(as.name(paste("dataset_subset",i, sep=""))), ti_method())
    temp_model <- infer_trajectory(add_prior_information(eval(as.name(paste("dataset_subset",i, sep=""))), start_id ='ACTATCTGTGTTTGGT-1'), ti_method())
    temp_dr_model <- temp_model %>% add_dimred(dyndimred::dimred_pca, expression_source = eval(as.name(paste("dataset_subset",i, sep=""))))
    temp_dr_model <- add_cell_waypoints(temp_dr_model)
    do.call("<-",list(ti_dr_model, temp_dr_model))
  }
}

method <- "wishbone"
# Calculate three metrics and aggregate score
wishbone_results_2 <- matrix(nrow = 1, ncol = 3)
for (i in 1:3){
  # Calculate correlation
  correlation <- calculate_metrics(eval(as.name(paste(method,"_dr_model", sep=""))),eval(as.name(paste(method,"_dr_model_subset",i, sep=""))), metrics$metric_id[1])$correlation
  # Calculate normalized rf_mse
  rf_mse <- calculate_metrics(eval(as.name(paste(method,"_dr_model", sep=""))),eval(as.name(paste(method,"_dr_model_subset",i, sep=""))), metrics$metric_id[3])$rf_nmse
  # Calculate edge-flip score
  edge_flip <- calculate_metrics(eval(as.name(paste(method,"_dr_model", sep=""))),eval(as.name(paste(method,"_dr_model_subset",i, sep=""))), metrics$metric_id[8])$edge_flip
  # Record harmonic mean of scores
  wishbone_results_2[1,i] <- calculate_harmonic_mean(correlation,rf_mse,edge_flip)
}
# Calculate total score (arithmetic mean)
wishbone_score[1,2] = mean(wishbone_results_2)

# Save as png with 1300/720 aspect ratio
plot <- "phenopath"
patchwork::wrap_plots(
  plot_dimred(eval(as.name(paste(plot,"_dr_model", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_base", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset1", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_1", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset2", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_2", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset3", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_3", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset4", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_4", sep=""))
)

# Re-dimred for better graphics
#paga_tree_dr_model_subset1 <- paga_tree_dr_model_subset1 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset1$expression)
#paga_tree_dr_model_subset2 <- paga_tree_dr_model_subset2 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset2$expression)
#paga_tree_dr_model_subset3 <- paga_tree_dr_model_subset3 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset3$expression)
#phenopath_dr_model_subset4 <- phenopath_dr_model_subset4 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset4$expression)
