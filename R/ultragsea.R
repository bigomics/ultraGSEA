##---------------------------------------------------------------
##-------------------- ultragsea --------------------------------
##---------------------------------------------------------------

#' Ultra-fast GSEA using z-score or geneset correlation. Results of
#' these methods highly correlate with GSEA/fGSEA but are much faster.
#'
#' @export
ultragsea <- function(G, fc, alpha=0, minLE=1, corshrink=3,
                      minsize = 1L, maxsize = 9999L, center=TRUE,
                      method=c("cor","ztest","ttest","goat","camera")[1],
                      format=c("simple","as.gsea")[1]) {

  gg <- intersect(names(fc), rownames(G))
  fc <- fc[gg]
  G <- G[gg,]
  size <- Matrix::colSums(G!=0)
  if(method == "goat") minsize <- max(minsize,10L)
  sel <- which(size >= minsize & size <= maxsize)
  G <- G[, sel, drop=FALSE]
  size <- size[sel]

  gmt <- NULL
  if(method %in% c("goat","camera") && is.null(gmt)) {
    gmt <- mat2gmt(G)
  }
  
  addLE <- (format=="as.gsea")
  p_value <- NULL
  q_value <- NULL
  stat_value <- NULL
  if(method == "ztest") {
    zres <- gset.ztest(G, fc, alpha=alpha, center=center)
    p_value <- zres$p
    stat_value <- zres$stat
  } else if(method == "ttest") {
    zres <- gset.ztest(G, fc, alpha=alpha, pdist="norm", center=center)
    p_value <- zres$p
    stat_value <- zres$stat
  } else if(method == "cor") {
    res <- gset.cortest(G, fc, compute.p=TRUE, use.rank=FALSE,
      corshrink = corshrink) 
    p_value <- res$p.value[,1]
    q_value <- res$q.value[,1]
    stat_value <- res$rho[,1]
  } else if(method=="goat") {
    if(!require("goat")) stop("please install goat")    
    res <- goat(pathways=gmt, stats=fc, filter=FALSE, method="goat") 
    res <- res[match(colnames(G),res$pathway),]
    p_value <- res$pval
    q_value <- res$padj
    stat_value <- res$score
  } else if(method=="camera") {
    if(!require("limma")) stop("please install limma")
    res <- limma::cameraPR(fc, gmt)
    res <- res[match(colnames(G),rownames(res)),]
    p_value <- res$PValue
    q_value <- res$FDR
    stat_value <- -log10(res$PValue)*(1 - 2*(res$Direction=="Down"))
  } else {
    stop("unsupported method: ", method)
  }

  ## Add leading edge list. The leading edge list is computed as those
  ## genes in the gene set that have fold-change larger than z score
  ## z>minLE (default minLE=1).
  leading_edge <- NULL
  if(addLE) {
    zsign <- sign(stat_value) 
    zmat <- sweep((G!=0), 2, zsign, FUN='*')
    zmat <- (zmat * fc > minLE) ## signed            
    leading_idx <- Matrix::which(zmat == 1, arr.ind=TRUE)  
    leading_edge <- tapply( leading_idx[,1], leading_idx[,2], list)
    leading_edge <- lapply( leading_edge, function(i) rownames(zmat)[i])
    names(leading_edge) <- colnames(zmat)[as.integer(names(leading_edge))]
    leading_edge <- leading_edge[match(colnames(zmat),names(leading_edge))]
    names(leading_edge) <- colnames(zmat)
  }
  
  #size <- Matrix::colSums(G!=0)
  if(is.null(q_value)) {
    q_value <- p.adjust(p_value, method="fdr")
  }

  if(format == "as.gsea") {
    # normalize like GSEA
    stat_value0 <-  stat_value / max(abs(stat_value),na.rm=TRUE) 
    stat_value1 <-  stat_value / sd(stat_value,na.rm=TRUE)
    df <- data.table::data.table(
      pathway = colnames(G),
      pval = p_value,
      padj = q_value,
      log2err = as.numeric(NA),
      ES = stat_value0,
      NES = stat_value1,
      size = size,
      leadingEdge = leading_edge        
    )
  } else {
    ## simple format
    df <- data.frame(
      pathway = colnames(G),
      pval = p_value,
      padj = q_value,
      score = stat_value,
      size = size
    )
    rownames(df) <- names(p_value)
  }
  df
}

#' Fast replacement of fgsea with p/q values computed with ultragsea.cor
#' 
#' @export
fgsea <- function(pathways, stats, minSize = 1, maxSize = length(stats)-1,
                  nperm=20) {
  G <- gmt2mat(pathways)
  gg <- intersect(names(stats),rownames(G))
  G <- G[gg,match(names(pathways),colnames(G))]
  stats <- stats[gg]
  gsize <- Matrix::colSums(G)
  G <- G[, gsize >= minSize & gsize <= maxSize]
  # run fgseaSimple only for few iterations only for NES
  suppressWarnings(fres <- fgsea::fgseaSimple(pathways, stats, nperm=nperm))
  # use ultragsea for p-value
  res <- ultragsea(G, stats, method="cor", alpha=0, format="simple")
  res <- res[match(res$pathway,fres$pathway),]  
  fres$nMoreExtreme <- NULL
  fres$pval <-  res$pval
  fres$padj <-  res$padj  
  fres
}

