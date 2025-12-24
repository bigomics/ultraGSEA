
run_plaid = function(se, matG){
  
  es = plaid::plaid(assay(se), matG)

  # set design matrix
  group = factor(se$GROUP)
  design = model.matrix(formula("~group"))  
  
  # fit the linear model to the GSVA enrichment scores
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
