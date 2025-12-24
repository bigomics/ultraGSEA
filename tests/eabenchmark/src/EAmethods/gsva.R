
run_gsva = function(se, kegg_list, matG, seed_n, use.plaid=FALSE){
 
  if(use.plaid) {
    es <- plaid::replaid.gsva( assay(se), matG)
  } else {
    kcdf = ifelse(se@metadata$dataType=='rseq', "Poisson", "Gaussian")
    set.seed(seed_n)
    es = GSVA::gsva(expr=assay(se), 
      gset.idx.list=kegg_list, 
      kcdf=kcdf,
      parallel.sz=2)
  }
  
  # fit the linear model to the GSVA enrichment scores
  group = factor(se$GROUP)
  design = model.matrix(formula("~group"))  
  fit = limma::lmFit(es, design)
  fit = limma::eBayes(fit)
  res = limma::topTable(fit, 
                        number=nrow(es), 
                        coef="group1", 
                        sort.by="none", 
                        adjust.method="none")
  
  # process output
  res = res[,c("t", "P.Value")]
  return(res)
}
