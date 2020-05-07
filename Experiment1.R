guidelines_shiny(dataset)

# For each fo the methods selected, create a TI model, 
# then apply the same dimensionality reduction (pca) to each of the models

# run with add_prior_information(dataset, start_id ='ACTATCTGTGTTTGGT-1')
prior_methods <- c("paga_tree","paga","wanderlust","wishbone")

methods <- c("slingshot","scorpius","angle","embeddr","mst","waterfall","tscan","comp1","slice","elpilinear",
             "phenopath","pcreode","elpicycle","celltree_maptpx","dpt","celltree_vem","ouijaflow")
for (method in methods){
  for (i in 1:4){
    # Model and function variable names
    ti_dr_model <- paste(method,"_dr_model_",i, sep="")
    ti_method <- match.fun(paste("ti_",method, sep=""))
    # Run method
    #temp_model <- infer_trajectory(add_prior_information(dataset, start_id ='ACTATCTGTGTTTGGT-1'), ti_method())
    temp_model <- infer_trajectory(dataset, ti_method())
    temp_dr_model <- temp_model %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset)
    temp_dr_model <- add_cell_waypoints(temp_dr_model)
    do.call("<-",list(ti_dr_model, temp_dr_model))
  }
}

# Calculate three metrics and aggregate score
model <- "wishbone"
wishbone_results_1 <- matrix(nrow = 5, ncol = 5)
wishbone_score <- matrix(nrow=1, ncol=3)
wishbone_dr_model_5 <- wishbone_dr_model
for (i in 1:5){
  for (j in 1:5){
    # A model has a similarity score of 1 with itself
    if (i == j){
      wishbone_results_1[i,j] <- 1
    } 
    # We have a symmetric matrix, so just copy the values over
    else if (i > j){
      wishbone_results_1[i,j] <- wishbone_results_1[j,i]
    } 
    else {
      # Calculate correlation
      correlation <- calculate_metrics(eval(as.name(paste(model,"_dr_model_",i, sep=""))),eval(as.name(paste(model,"_dr_model_",j, sep=""))), metrics$metric_id[1])$correlation
      # Calculate normalized rf_mse
      rf_mse <- calculate_metrics(eval(as.name(paste(model,"_dr_model_",i, sep=""))),eval(as.name(paste(model,"_dr_model_",j, sep=""))), metrics$metric_id[2])$rf_mse
      rf_mse <- (rf_mse-0.3)/-0.3
      # Calculate edge-flip score
      edge_flip <- calculate_metrics(eval(as.name(paste(model,"_dr_model_",i, sep=""))),eval(as.name(paste(model,"_dr_model_",j, sep=""))), metrics$metric_id[8])$edge_flip
      # Record harmonic mean of scores
      wishbone_results_1[i,j] <- calculate_harmonic_mean(correlation,rf_mse,edge_flip)
    }
  }
}
# Calculate total score (arithmetic mean excluding the diagonal)
wishbone_score[1,1] = mean(wishbone_results_1[row(wishbone_results_1)!=col(wishbone_results_1)])

# Save as png with 1300/720 aspect ratio
plot <- "celltree_maptpx"
patchwork::wrap_plots(
  plot_dimred(eval(as.name(paste(plot,"_dr_model", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_0", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_1", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_1", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_2", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_2", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_3", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_3", sep="")),
  plot_dimred(eval(as.name(paste(plot,"_dr_model_4", sep=""))), color_cells = "pseudotime") + ggtitle(paste(plot,"_4", sep=""))
)

# Test which TI can be run
#slicer_model <- infer_trajectory(add_prior_information(dataset, start_id ='ACTATCTGTGTTTGGT-1') , ti_slicer())
#slicer_dr_model <- slicer_model %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset$expression)
#slicer_dr_model <- add_cell_waypoints(slicer_dr_model)
#plot_dimred(slicer_dr_model, color_cells = "pseudotime") + ggtitle("SLICER")

# Re-dimred for better graphics
#celltree_maptpx_dr_model_1 <- celltree_maptpx_dr_model_1 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset$expression)
#celltree_maptpx_dr_model_2 <- celltree_maptpx_dr_model_2 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset$expression)
#celltree_maptpx_dr_model_3 <- celltree_maptpx_dr_model_3 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset$expression)
#celltree_maptpx_dr_model_4 <- celltree_maptpx_dr_model_4 %>% add_dimred(dyndimred::dimred_pca, expression_source = dataset$expression)


