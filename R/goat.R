goat_test <- function(fc, gmt, method=c("goat","gsea")[1]) {

  genesets  <- data.table::data.table(
    source = "GSET", source_version = "local",
    id = paste0("gset",1:length(gmt)),
    name = names(gmt),
    genes = gmt,
    ngenes = sapply(gmt,length),
    ngenes_signif = 0L
  )
  genesets <- tibble::as_tibble(genesets)
  
  genelist <- tibble::tibble(
    gene = names(fc),
    test = FALSE,
    signif = FALSE,  
    effectsize = fc
  )
  
  result <- goat::test_genesets(
    genesets, genelist,
    method = method,
    score_type = "effectsize",
    padj_method = "BH",
    padj_cutoff = 0.05
  )
  if(method == "gsea") score <- result$gsea_nes
  if(method == "goat") score <- result$zscore
  
  df <- data.frame(
    pathway = result$name,
    pval = result$pvalue,
    padj = result$pvalue_adjust,
    score = score,  
    size = result$ngenes
  )
  df
}

