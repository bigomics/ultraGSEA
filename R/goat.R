goat <- function(pathways, stats, minSize = 10L, maxSize = 1500L,
                 filter = TRUE, method=c("goat","gsea")[1]) {
  method <- method[1]

  if(method %in% c("fisherexact","fisherexact_ease","hypergeometric")) {
    stop("method not supported. please use goat::test_genesets(...)")
  }
  if(!method %in% c("goat","goat_fitfunction","goat_precomputed",
    "goat_bootstrap","gsea")) {
    stop("invalid method")
  }
  
  ## goat does not like length(stats) larger than 20000
  if(length(stats) > 20000 && method %in% c("goat","goat_precomputed")) {
    message("truncating stats vector")
    sel <- head(order(-abs(stats)),20000)
    stats <- stats[sel]
    filter <- TRUE ## need filter
  }
  
  genesets  <- data.table::data.table(
    source = "SOURCE",
    source_version = paste("goat",goat::goat_version()),
    id = paste0("pathway",1:length(pathways)),
    name = names(pathways),
    genes = pathways,
    ngenes = sapply(pathways,length),
    ngenes_signif = 0L
  )
  genesets <- tibble::as_tibble(genesets)
  
  genelist <- tibble::tibble(
    gene = names(stats),
    test = FALSE,
    signif = FALSE,  
    effectsize = stats
  )

  if(filter) {
    message("filtering genesets...")
    ## filter geneset sizes. For speed, disable if already filtered
    genesets <- goat::filter_genesets(
      genesets,
      genelist,
      min_overlap = minSize,
      max_overlap = maxSize
    )
  }
  
  result <- goat::test_genesets(
    genesets, genelist,
    method = method,
    score_type = "effectsize",
    padj_method = "BH",
    padj_cutoff = 0.05
  )

  if(method == "gsea") {
    score <- result$gsea_nes
  } else {
    score <- result$zscore
  }
  
  df <- data.frame(
    pathway = result$name,
    pval = result$pvalue,
    padj = result$pvalue_adjust,
    score = score,  
    size = result$ngenes
  )
  df
}

