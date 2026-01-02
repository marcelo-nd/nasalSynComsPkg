utils::globalVariables(c(
  # dplyr / tidy evaluation
  ".data", ".id",
  "sample_id", "syncom_id", "time_label", "time_num", "rep_raw",
  "syncom_order", "pairwise_dist", "dist_to_TF", "mean_dist_to_TF",
  "time_raw",
  
  # plot columns / aesthetics seen in your functions
  "abundance", "species", "species2", "strain", "experiment",
  "Bacteria", "Sample", "Abundance", "Metabolite", "Replicate",
  "Relative", "log2FC", "ctrl_mean", "TargetCluster", "adj.P.Val", "logFC",
  "PCo1", "PCo2", "taxonomy"
))