rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='1',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='fgsea',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores

#################################################
## FUNTIONAL GENE SETS
kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff')
kegg_df = kegg_df %>%
  filter(pathway_id %in%
           as.character(
             kegg_df %>%
               select(pathway_id,entrezgene_id) %>%
               unique() %>%
               group_by(pathway_id) %>%
               summarise(n=n()) %>%
               filter(n>14) %>%
               droplevels() %>%
               pull(pathway_id))) %>%
  droplevels() %>%
  select(pathway_id,entrezgene_id) %>%
  dplyr::rename(GeneID=entrezgene_id)

kegg_list = tapply(kegg_df$GeneID, kegg_df$pathway_id, function(x) as.character(x))
kegg_list = lapply(kegg_list, function(x) unique(x))
kegg_list = kegg_list[sapply(kegg_list,length)>=10]

matG <- ultragsea::gmt2mat(kegg_list)

#################################################
## EXPRESSION DATA
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
                      col_types = 'fcfffiddiciccccc') %>% 
  dplyr::rename(seed_n = seed.number)
head(dataset_df)

outdir = paste0('results/TPbenchmark/')

#################################################
## EXECUTE METHOD
run_enrichment = function(d){
  library(tidyverse)
  message(d)
  #################################################
  ## EXPRESSION DATA
  #x = dataset_df[1,] %>% droplevels()
  x = dataset_df[dataset_df$GEO==d,] %>%
    droplevels()
  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
  se@colData$GROUP = as.numeric(as.character(se@colData$GROUP))
  se@metadata$dataType = as.character(x[['DataType']])
  print(se)
  if (method_name == 'camera'){ 
   ##########################################
   #### CAMERA
   source('src/EAmethods/camera.R')
   ptm = proc.time()
   tmp_enrichment_df = 
     tryCatch(
       run_camera(se, kegg_list, var_corr = FALSE),
       error=function(cond) {
         return(NULL)})
   if (!is.null(tmp_enrichment_df)){
     p = rownames(tmp_enrichment_df)
     tmp_enrichment_df = tmp_enrichment_df %>%
       dplyr::rename(Pvalue=PValue) %>%
       magrittr::set_rownames(NULL) %>%
       add_column(Correlation=0.01) %>%
       add_column(pathway_id = p,
                  GEO = x[['GEO']]) %>%
       select(GEO,pathway_id,NGenes,Direction,Pvalue,Correlation)
   } else {return()}
   
  } else if (method_name == 'camera*'){
    ##########################################
    #### CAMERA*
    source('src/EAmethods/camera.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_camera(se, kegg_list, var_corr = TRUE),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,NGenes,Direction,Pvalue,Correlation)
    } else {return()}
    
  } else if (method_name == 'fgsea'){
    ##########################################
    #### fGSEA
    source('src/EAmethods/fgsea.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_fgseaMultilevel(se, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,ES,NES,Pvalue,size)
    } else {return()}

  } else if (grepl("^ultragsea",method_name)){
    ##########################################
    #### ultraGSEA
    source('src/EAmethods/ultragsea.R')
    ptm = proc.time()
    ## method_name='ultragsea.cor.a0.c1'
    pars <- strsplit(method_name,split="[.]")[[1]]
    method="cor";alpha=0.5;center=1;use.rank=FALSE
    if(length(pars)>=2) method <- pars[2]
    if(length(pars)>=3) alpha  <- as.numeric(sub("^a","",pars[3]))
    if(length(pars)>=4) center <- as.integer(sub("^c","",pars[4]))
    if(length(pars)>=5 && pars[5]=="r") use.rank <- TRUE

    tmp_enrichment_df = tryCatch(
      run_ultragsea(se, matG, method=method, alpha=alpha,
        center=center, use.rank=use.rank),
      error=function(cond) {
        return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,score,Pvalue,size)
    } else {return()}
    
  } else if (method_name == 'goat'){
    ##########################################
    #### ultraGSEA
    source('src/EAmethods/goat.R')
    ptm = proc.time()
    tmp_enrichment_df = tryCatch(
      run_goat(se, kegg_list),
      error=function(cond) {
        return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,Pvalue)
    } else {return()}
  } else if (method_name == 'cameraPR'){
    ##########################################
    #### cameraPR
    source('src/EAmethods/cameraPR.R')
    ptm = proc.time()
    tmp_enrichment_df = tryCatch(
      run_cameraPR(se, kegg_list),
      error=function(cond) {
        return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = PValue) %>%
        add_column(GEO = x[['GEO']]) %>%
        add_column(pathway_id = names(kegg_list)) %>%
        select(GEO,pathway_id,Pvalue)
    } else {return()}
    
  } else if (method_name == 'gsea'){
    ##########################################
    #### GSEA
    source('src/EAmethods/gsea.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_gsea(se, kegg_df, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df@result %>%
        dplyr::rename(pathway_id = ID,
                      Pvalue = pvalue) %>%
        add_column(GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,enrichmentScore,NES,Pvalue,setSize)
    } else {return()}
    
  } else if (method_name == 'gsa'){
    ##########################################
    #### GSA
    source('src/EAmethods/gsa.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
         run_gsa(se, kegg_list, x[['seed_n']]),
         error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        add_column(pathway_id = names(kegg_list),
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,ES,Pvalue)
    } else {return()}
    
  } else if (method_name %in% c('gsva','replaid.gsva')){
    ##########################################
    #### GSVA
    source('src/EAmethods/gsva.R')
    ptm = proc.time()
    use.plaid <- (method_name == "replaid.gsva")
    tmp_enrichment_df = 
      tryCatch(
        run_gsva(se, kegg_list, matG, x[['seed_n']],
          use.plaid = use.plaid),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = P.Value) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,t,Pvalue)
    } else {return()}

  } else if (method_name == 'plaid'){
    ##########################################
    #### PLAID
    source('src/EAmethods/plaid.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_plaid(se, matG),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = P.Value) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,t,Pvalue)
    } else {return()}
    
  } else if (method_name == 'padog'){
    ##########################################
    #### PADOG
    source('src/EAmethods/padog.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_padog(se, kegg_list, 8, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      # p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        magrittr::set_rownames(NULL) %>%
        dplyr::rename(Pvalue = Ppadog,
                      pathway_id = ID) %>%
        add_column(GEO = x[['GEO']]) %>%
        mutate(meanAbsT0 = round(meanAbsT0,2),
               padog0 = round(padog0,4)) %>%
        select(GEO,pathway_id,meanAbsT0,padog0,PmeanAbsT,Pvalue)
    } else {return()}
    
  } else if (method_name == 'roast'){
    ##########################################
    #### ROAST
    source('src/EAmethods/roast.R')
    ptm = proc.time()
    tmp_enrichment_df = 
      tryCatch(
        run_roast(se, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = x[['GEO']]) %>%
        select(GEO,pathway_id,PropDown,PropUp,Direction,Pvalue)
    } else {return()}
  }
  
  runtime = proc.time() - ptm
  tmp_enrichment_df$user=round(runtime[[1]],2)
  tmp_enrichment_df$system=round(runtime[[2]],2)
  tmp_enrichment_df$elapsed=round(runtime[[3]],2)
  ### FOR ALL METHODS
  return(tmp_enrichment_df)
}

do.parallel=FALSE
if(do.parallel) {
  pb = utils::txtProgressBar(min=0, max=nrow(dataset_df), style = 3)
  progress = function(n) utils::setTxtProgressBar(pb, n)
  opts = list(progress = progress)
  cl = parallel::makeCluster(ncores)
  doSNOW::registerDoSNOW(cl)
  boot = foreach::foreach(i = dataset_df$GEO, .options.snow = opts)
  output_df = foreach::`%dopar%`(boot, run_enrichment(i)) %>%
    data.table::rbindlist()
  parallel::stopCluster(cl)
} else {
  output_df <- lapply(dataset_df$GEO, function(i) run_enrichment(i))
  output_df <- output_df %>% data.table::rbindlist()
}

## output TPbenchmark
output_df %>%
  select(-c(user,system,elapsed)) %>%
  write_tsv(file = paste0(outdir,method_name,'.tsv.gz'))

## output runtime
output_df %>%
  select(GEO,user,system,elapsed) %>%
  unique() %>%
  write_tsv(file = paste0('results/runtime/',method_name,'.tsv.gz'))
