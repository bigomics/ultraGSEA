
run_ultragsea = function(se, matG, method, alpha, center,
                         use.rank = FALSE){
  ## if nperm!=NULL, run fgseaSimple, otherwise run fgseaMultilevel
  dea_df = as.data.frame(SummarizedExperiment::rowData(se)) 
  dea_df = dea_df[order(dea_df$t, decreasing = TRUE),]
  dea_ranks = dea_df$t
  names(dea_ranks) = rownames(dea_df)
  dea_ranks = dea_ranks[!is.na(dea_ranks)]
  if(use.rank) dea_ranks <- rank(dea_ranks)
  ultragsea::ultragsea(fc = dea_ranks, G = matG, method = method,
    alpha = alpha, center = center)
}

