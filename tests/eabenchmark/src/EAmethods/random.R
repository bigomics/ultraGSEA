#!/usr/bin/env Rscript

run_random = function(genes,gene_sets_list){                          
  
  ## Extract pathway ids, number of pahtways
  gs_names = names(gene_sets_list)
  gs_length = length(gene_sets_list)
  
  ## Sample pvalues from uniform distribution
  set.seed(length(genes))
  pval_uniform = runif(gs_length,min=0,max=1)
  
  ## Output a data frame with pathway_id,Pvalue as columns
  enrichment_df = data.frame(
    pathway_id = gs_names,
    Pvalue = pval_uniform
    )  
   
  return(enrichment_df)
}
