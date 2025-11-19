#' Utility to convert fgsea/ultragsea style results to 'gseaResult'
#' object compatible with clusterProfiler and enrichplot
#' functions. Using this function you can use all plotting functions
#' from enrichplot.
#'
#' @export
new.gseaResult <- function(res, geneList, geneSets) 
{
  req.cols <- c("pathway","NES","pval","leadingEdge")
  if(!all(req.cols %in% colnames(res))) stop("missing columns in result table")
  params <- list(pvalueCutoff = 0.05, nPerm = 1000, 
    pAdjustMethod = "fdr", exponent = 1, minGSSize = 1, maxGSSize = 99999)
  ##geneSets <- lapply(geneSets, function(x) intersect(x,names(geneList)))
  ID = as.character(res$pathway)
  if(is.null(res$size)) {
    res$size <- sapply(geneSets,length)[match(ID,names(geneSets))]
  }
  if(is.null(res$leadingEdge)) {
    sdfc <- 1.0 * sd(fc, na.rm=TRUE)
    sig.genes <- names(geneList[abs(geneList) > sdfc])
  }
  if(is.null(res$ES)) res$NES
  if(is.null(res$padj)) p.adjust(res$pval, method="fdr")
  df <- data.frame(
    ID = as.character(res$pathway),
    Description = res$pathway, 
    setSize = res$size, enrichmentScore = res$ES, 
    NES = res$NES, pvalue = res$pval, p.adjust = res$padj, 
    qvalue = res$padj, stringsAsFactors = FALSE)
  row.names(df) <- df$ID
  numle <- sapply(res$leadingEdge,length)
  tags <- paste0(round(100*(numle/df$setSize)),"%")
  df$rank <- NA
  df$leading_edge <- tags
  df$core_enrichment <- sapply(res$leadingEdge,paste,collapse = "/")
  geneList <- sort(geneList, decreasing=TRUE)
  new("gseaResult", result = df, geneSets = geneSets, geneList = geneList, 
    params = params, readable = FALSE)
}
