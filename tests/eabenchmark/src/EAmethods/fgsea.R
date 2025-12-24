
run_fgseaSimple = function(se,
                    kegg_list,
                    seed_n){
  ## if nperm!=NULL, run fgseaSimple, otherwise run fgseaMultilevel
  dea_df = as.data.frame(SummarizedExperiment::rowData(se)) 
  dea_df = dea_df[order(dea_df$t, decreasing = T),]
  dea_ranks = dea_df$t
  names(dea_ranks) = rownames(dea_df)
  dea_ranks = dea_ranks[!is.na(dea_ranks)]
  
  set.seed(seed_n)
  fgsea::fgsea(pathways = kegg_list, 
               stats = dea_ranks,
               nperm = 2000,
               minSize = 1,
               maxSize = 1000)
}

run_fgseaMultilevel = function(se,
                               kegg_list,
                               seed_n){
  ## if nperm!=NULL, run fgseaSimple, otherwise run fgseaMultilevel
  dea_df = as.data.frame(SummarizedExperiment::rowData(se)) 
  dea_df = dea_df[order(dea_df$t, decreasing = T),]
  dea_ranks = dea_df$t
  names(dea_ranks) = rownames(dea_df)
  dea_ranks = dea_ranks[!is.na(dea_ranks)]
  
  set.seed(seed_n)
  fgsea::fgsea(pathways = kegg_list, 
               stats = dea_ranks,
               eps = 0.0,
               minSize = 1,
               maxSize = 1000,
               nPermSimple = 2000)
}