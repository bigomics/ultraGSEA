
run_cameraPR = function(se, kegg_list){
  dea_df = as.data.frame(SummarizedExperiment::rowData(se)) 
  dea_df = dea_df[order(dea_df$t, decreasing = TRUE),]
  dea_ranks = dea_df$t
  names(dea_ranks) = rownames(dea_df)
  dea_ranks = dea_ranks[!is.na(dea_ranks)]
  limma::cameraPR(dea_ranks, kegg_list)
}

