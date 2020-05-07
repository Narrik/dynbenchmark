# Set seed and randomly shuffle matrix
set.seed(42)
# Set prop=0.125 for extra genes
hvg_extra <- getTopHVGs(dec, prop=0.125)
sce_counts_extra <- counts(sce)[hvg_extra,]
sce_logcounts_extra <- logcounts(sce)[hvg_extra,]
cols <- sample(ncol(countsMatrix))

# Transform from SCE to matrix
countsMatrix_extra <- sce_counts_extra %>% as.matrix() %>% t
logCountsMatrix_extra <- sce_logcounts_extra %>% as.matrix() %>% t

countsMatrix_subset1_extra <- countsMatrix_extra[,c(cols[373:1485],1486:1857)]
countsMatrix_subset2_extra <- countsMatrix_extra[,c(cols[c(1:372,745:1485)],1486:1857)]
countsMatrix_subset3_extra <- countsMatrix_extra[,c(cols[c(1:744,1117:1485)],1486:1857)]
countsMatrix_subset4_extra <- countsMatrix_extra[,c(cols[1:1113],1486:1857)]

logCountsMatrix_subset1_extra <- logCountsMatrix_extra[,c(cols[373:1485],1486:1857)]
logCountsMatrix_subset2_extra <- logCountsMatrix_extra[,c(cols[c(1:372,745:1485)],1486:1857)]
logCountsMatrix_subset3_extra <- logCountsMatrix_extra[,c(cols[c(1:744,1117:1485)],1486:1857)]
logCountsMatrix_subset4_extra <- logCountsMatrix_extra[,c(cols[1:1113],1486:1857)]

dataset_subset1_extra<- wrap_expression(counts = countsMatrix_subset1_extra,expression = logCountsMatrix_subset1_extra)
dataset_subset2_extra<- wrap_expression(counts = countsMatrix_subset2_extra,expression = logCountsMatrix_subset2_extra)
dataset_subset3_extra<- wrap_expression(counts = countsMatrix_subset3_extra,expression = logCountsMatrix_subset3_extra)
dataset_subset4_extra<- wrap_expression(counts = countsMatrix_subset4_extra,expression = logCountsMatrix_subset4_extra)

methods <- c("comp1")
for (method in methods){
  for (i in 1:4){
    # Model and function variable names
    ti_dr_model <- paste(method,"_dr_model_subset",i,"_extra", sep="")
    ti_method <- match.fun(paste("ti_",method, sep=""))
    # Run method
    temp_model <- infer_trajectory(eval(as.name(paste("dataset_subset",i,"_extra", sep=""))), ti_method())
    #temp_model <- infer_trajectory(add_prior_information(eval(as.name(paste("dataset_subset",i,"_extra", sep=""))), start_id ='ACTATCTGTGTTTGGT-1'), ti_method())
    temp_dr_model <- temp_model %>% add_dimred(dyndimred::dimred_pca, expression_source = eval(as.name(paste("dataset_subset",i,"_extra", sep=""))))
    temp_dr_model <- add_cell_waypoints(temp_dr_model)
    do.call("<-",list(ti_dr_model, temp_dr_model))
  }
}

2# Calculate three metrics and aggregate score
wishbone_results_3 <- matrix(nrow = 1, ncol = 4)
for (i in 1:4){
  # Calculate correlation
  correlation <- calculate_metrics(eval(as.name(paste(method,"_dr_model", sep=""))),eval(as.name(paste(method,"_dr_model_subset",i,"_extra", sep=""))), metrics$metric_id[1])$correlation
  # Calculate anormalized rf_mse
  rf_mse <- calculate_metrics(eval(as.name(paste(method,"_dr_model", sep=""))),eval(as.name(paste(method,"_dr_model_subset",i,"_extra", sep=""))), metrics$metric_id[3])$rf_nmse
  # Calculate edge-flip score
  edge_flip <- calculate_metrics(eval(as.name(paste(method,"_dr_model", sep=""))),eval(as.name(paste(method,"_dr_model_subset",i,"_extra", sep=""))), metrics$metric_id[8])$edge_flip
  # Record harmonic mean of scores
  wishbone_results_3[1,i] <- calculate_harmonic_mean(correlation,rf_mse,edge_flip)
}
# Calculate total score (arithmetic mean)
wishbone_score[1,3] = mean(wishbone_results_3)

# Save as png with 1300/720 aspect ratio
plot <- "comp1"
patchwork::wrap_plots(
  plot_dimred(eval(as.name(paste(plot,"_dr_model", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_base", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset1_extra", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_1", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset2_extra", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_2", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset3_extra", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_3", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_subset4_extra", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_4", sep=""))
)

# Re-dimred for better graphics
#paga_dr_model_subset1_extra <- paga_dr_model_subset1_extra %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset1_extra$expression)
#paga_dr_model_subset2_extra <- paga_dr_model_subset2_extra %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset2_extra$expression)
#comp1_dr_model_subset3_extra <- comp1_dr_model_subset3_extra %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset3_extra$expression)
#comp1_dr_model_subset4_extra <- comp1_dr_model_subset4_extra %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset_subset4_extra$expression)
