#!/usr/bin/env Rscript

run_fastFET = function(genes, matG, bg){                          
  res <- ultragsea::gset.fastFET(genes, matG, bg, report.genes=FALSE) 
  return(res)
}

run_fisher = function(A,PT=18000,gs){                          
  # LH: list hits
  # LT: list total
  # PH: population hits
  # PT: population total
  # User genes | genome
  # ------------------
  # LH    | PH-LH | PH
  # ------------------
  # LT-LH   | PT-LT-(PH-LH) | PT-PH
  # ------------------
  
  A = as.character(A)
  LT = length(A)
  PT = as.numeric(PT) 
  
  gs_names = names(gs)
  enrichment_df = lapply(gs_names, function(gs_n){
    B = as.character(gs[[gs_n]])
    PH = length(B)
    LH = length(intersect(A,B))

    contingency_table = rbind(c(LH, PH-LH),
                              c(LT-LH,PT-LT-(PH-LH))) 
    
    pval=fisher.test(contingency_table, alternative = 'greater')[['p.value']]
    
    return(
      data.frame(
        pathway_id = gs_n,
        Overlap = LH,
        Pvalue = pval,
        GeneID = ifelse(LH>0,paste0(intersect(A,B), collapse = ','),NA))
    )   
  }) %>% data.table::rbindlist()
  return(enrichment_df)
}
