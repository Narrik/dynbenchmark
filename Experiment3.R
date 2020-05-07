# Calculate three metrics and aggregate score
methods <- c("slingshot","paga_tree","scorpius","angle","paga","embeddr","mst","waterfall","tscan","comp1","slice",
             "elpilinear","phenopath","wanderlust","wishbone","elpicycle","celltree_maptpx","dpt","celltree_vem","ouijaflow")
final_results <- matrix(nrow = 20, ncol = 20)
i <- 1
j <- 1
for (method1 in methods){
  for (method2 in methods){
    # A model has a similarity score of 1 with itself
    if (i == j){
      final_results[i,j] <- 1
    } 
    # We have a symmetric matrix, so just copy the values over
    else if (i > j){
      final_results[i,j] <- final_results[j,i]
    } 
    else {
      # Calculate correlation
      correlation <- calculate_metrics(eval(as.name(paste(method1,"_dr_model", sep=""))),eval(as.name(paste(method2,"_dr_model", sep=""))), metrics$metric_id[1])$correlation
      # Calculate normalized rf_mse
      rf_mse <- calculate_metrics(eval(as.name(paste(method1,"_dr_model", sep=""))),eval(as.name(paste(method2,"_dr_model", sep=""))), metrics$metric_id[3])$rf_nmse
      # Calculate edge-flip score
      edge_flip <- calculate_metrics(eval(as.name(paste(method1,"_dr_model", sep=""))),eval(as.name(paste(method2,"_dr_model", sep=""))), metrics$metric_id[8])$edge_flip
      # Record harmonic mean of scores
      final_results[i,j] <- calculate_harmonic_mean(correlation,rf_mse,edge_flip)
    }
    j <- j+1
  }
  i <- i+1
  j <- 1
}

# Save as png with 1300/720 aspect ratio
patchwork::wrap_plots(
  plot_dimred(scorpius_dr_model, color_cells = "pseudotime") + ggtitle("SCORPIUS"),
  plot_dimred(waterfall_dr_model, color_cells = "pseudotime") + ggtitle("Waterfall")
)
