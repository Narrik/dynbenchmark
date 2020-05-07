ranking <- c("slingshot","paga_tree","scorpius","angle","paga","embeddr","mst","waterfall","tscan","comp1","slice",
             "elpilinear","phenopath","wanderlust","wishbone","elpicycle","celltree_maptpx","dpt","celltree_vem","ouijaflow")
i <- 1
score1sum<- 0
score2sum<- 0
score3sum<- 0
score4sum<- 0
finalscoresum <- 0
for (method in ranking){
  score1 <- eval(as.name(paste(method,"_score",sep="")))[1,1]
  score1sum <- score1sum + score1
  score2 <- eval(as.name(paste(method,"_score",sep="")))[1,2]
  score2_scaled <- score2/0.6500408
  score2sum <- score2sum + score2_scaled
  score3 <- eval(as.name(paste(method,"_score",sep="")))[1,3]
  score3_scaled <- score3/0.9949668 
  score3sum <- score3sum + score3_scaled
  score4 <- (sum(final_results[,i])-1)/19
  score4_scaled <- score4/0.6748216
  score4sum <- score4sum + score4_scaled
  final_score <- calculate_harmonic_mean(score1,score2_scaled,score3,score4_scaled)
  finalscoresum <- finalscoresum + final_score
  #cat(method,score1,score2_scaled,score3_scaled,score4_scaled,final_score,"\n",sep=" & ")
  cat(method,round(score1,3),round(score2_scaled,3),round(score3_scaled,3),round(score4_scaled,3),round(final_score,3),"\n",sep=" & ")
  i <- i+1
}
 
# Save as png with 1300/720 aspect ratio
patchwork::wrap_plots(
  plot_dimred(slingshot_dr_model, color_cells = "pseudotime") + ggtitle("Slingshot"),
  plot_dimred(angle_dr_model, color_cells = "pseudotime") + ggtitle("Angle"),
  plot_dimred(embeddr_dr_model, color_cells = "pseudotime") + ggtitle("Embeddr"),
  plot_dimred(mst_dr_model, color_cells = "pseudotime") + ggtitle("MST"),
  plot_dimred(tscan_dr_model, color_cells = "pseudotime") + ggtitle("TSCAN"),
  plot_dimred(slice_dr_model, color_cells = "pseudotime") + ggtitle("SLICE")
)

# Save as png with 1300/720 aspect ratio
patchwork::wrap_plots(
  plot_dimred(elpilinear_dr_model, color_cells = "pseudotime") + ggtitle("ElPiGraph linear"),
  plot_dimred(wanderlust_dr_model, color_cells = "pseudotime") + ggtitle("Wanderlust"),
  plot_dimred(wishbone_dr_model, color_cells = "pseudotime") + ggtitle("Wishbone"),
  plot_dimred(elpicycle_dr_model, color_cells = "pseudotime") + ggtitle("ElPiGraph cycle"),
  plot_dimred(dpt_dr_model, color_cells = "pseudotime") + ggtitle("DPT"),
  plot_dimred(ouijaflow_dr_model, color_cells = "pseudotime") + ggtitle("ouijaflow")
)

sessionInfo()

