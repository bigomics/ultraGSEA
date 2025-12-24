
run_plaid = function(se, matG, method){
  
  X = assay(se)
  group = factor(se$GROUP)

  if(method == "ttest") {
    res <- ultragsea::plaid.ttest(X, matG, group, ref=0)
  }
  if(method == "limma") {
    res <- ultragsea::plaid.limma(X, matG, group, ref=0)
  }
  if(method == "cortest") {
    res <- ultragsea::plaid.cortest(X, matG, group, ref=0)
  }
  if(method == "dual") {
    res <- ultragsea::plaid.dualtest(X, matG, group, ref=0)
  }
    
  # process output
  res = res[,c("pval","padj")]
  return(res)
}
