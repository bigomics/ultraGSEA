rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='4',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-p", "--computer_name"), type="character", default='dbmini',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='padog',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores
message(paste0('Execute ', toupper(method_name), ' with #', ncores))

#################################################
## EXPRESSION DATA
dataset_df = read_tsv(file = paste0('data/FPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv'), col_types = 'fcfffiddiciccccc') %>%
  select(-c(GeneID,Ensembl,Ensembl_PRO)) %>%
  arrange(GEO)

dataset_df$GEO <- sub("[.].*","",dataset_df$GEO)
dataset_df <- dataset_df[!duplicated(dataset_df$GEO),]
dim(dataset_df)
d=dataset_df$GEO[1]

# %>% dplyr::rename(seed_n = seed.number)

outdir = 'results/FPbenchmark/'

#################################################
## EXECUTE METHOD
run_enrichment = function(d,iter,method_name){
  library(tidyverse)
  library(SummarizedExperiment)
  
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
  x = read_tsv(file = paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv'),
               col_types = 'fcfffiddiciccccc')
  d <- sub("[.].*","",as.character(d))
  x = read_tsv(file = paste0('data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv'),
               col_types = 'fcfffiddiciccccc') %>%
    select(-c(GeneID,Ensembl,Ensembl_PRO)) %>%
    filter(GEO==d) %>% 
    dplyr::rename(seed_n = seed.number) %>%
    droplevels()
  
  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
  se@colData$GROUP = as.numeric(as.character(se@colData$GROUP))
  se@metadata$dataType = as.character(x[['DataType']])
  
  permuted_se = as.character(
    read_tsv(file = paste0('data/FPbenchmark/',x[['GEO']],'.tsv.gz'), col_types = 'c', col_select = paste0('P',iter)) %>%
      pull(paste0('P',iter))
  )
  
  dat = assay(se)
  de = rowData(se)
  rownames(dat) = rownames(de) = permuted_se
  
  newSe =
    SummarizedExperiment(
      assays = as.matrix(dat),
      rowData = de,
      colData = colData(se),
      metadata = metadata(se)
    )
  newSe = subset(newSe,complete.cases(assay(newSe)))
  newSe@colData$GROUP = as.numeric(as.character(newSe@colData$GROUP))
  colnames(rowData(newSe)) = c('FC','p.value','t','ADJ.PVAL')
  newSe@metadata$dataType = x[['DataType']]
  
  rm(se)
  rm(de)
  rm(dat)
  rm(permuted_se)
  
  if (method_name == 'camera'){ 
    ##########################################
    #### CAMERA
    source('src/EAmethods/camera.R')
    tmp_enrichment_df = 
      tryCatch(
        run_camera(newSe, kegg_list, var_corr = FALSE),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(Correlation=0.01) %>%
        add_column(pathway_id = p,
                   GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,NGenes,Direction,Pvalue,Correlation)
    } else {return()}
    
  } else if (method_name == 'camera*'){
    ##########################################
    #### CAMERA*
    source('src/EAmethods/camera.R')
    tmp_enrichment_df = 
      tryCatch(
        run_camera(newSe, kegg_list, var_corr = TRUE),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,NGenes,Direction,Pvalue,Correlation)
    } else {return()}
    
  } else if (method_name == 'fgsea'){
    ##########################################
    #### fGSEA
    source('src/EAmethods/fgsea.R')
    tmp_enrichment_df = 
      tryCatch(
        run_fgseaMultilevel(newSe, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,ES,NES,Pvalue,size,log2err,leadingEdge)
    } else {return()}
    
  } else if (grepl("^ultragsea",method_name)){
    ##########################################
    #### ultraGSEA
    source('src/EAmethods/ultragsea.R')
    ## method_name='ultragsea.cor.a0.c1'
    pars <- strsplit(method_name,split="[.]")[[1]]
    method="cor";alpha=0;center=1;use.rank=FALSE
    if(length(pars)>=2) method <- pars[2]
    if(length(pars)>=3) alpha <- as.numeric(sub("^a","",pars[3]))
    if(length(pars)>=4) center <- as.integer(sub("^c","",pars[4]))
    if(length(pars)>=5 && pars[5]=="r") use.rank <- TRUE

    tmp_enrichment_df = tryCatch(
      run_ultragsea(newSe, matG, method=method, alpha=alpha,
        center=center, use.rank=use.rank),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,score,Pvalue,size)
    } else {return()}

  } else if (method_name == 'goat'){
    ##########################################
    #### ultraGSEA
    source('src/EAmethods/goat.R')
    tmp_enrichment_df = 
      tryCatch(
        run_goat(newSe, kegg_list),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(pathway_id = pathway,
                      Pvalue = pval) %>%
        add_column(GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,score,Pvalue,size)
    } else {return()}
    
  } else if (method_name == 'cameraPR'){
    ##########################################
    #### ultraGSEA
    source('src/EAmethods/cameraPR.R')
    tmp_enrichment_df = 
      tryCatch(
        run_cameraPR(newSe, kegg_list),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = PValue) %>%
        add_column(GEO = paste0(x[['GEO']],'.',iter)) %>%
        add_column(pathway_id = names(kegg_list)) %>%        
        select(GEO,pathway_id,Pvalue)
    } else {return()}
    
  } else if (method_name == 'gsea'){
    ##########################################
    #### GSEA
    source('src/EAmethods/gsea.R')
    tmp_enrichment_df = 
      tryCatch(
        run_gsea(newSe, kegg_df, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df@result %>%
        dplyr::rename(pathway_id = ID,
                      Pvalue = pvalue) %>%
        add_column(GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,enrichmentScore,NES,Pvalue,setSize,leading_edge,core_enrichment)
    } else {return()}
    
  } else if (method_name == 'gsa'){
    ##########################################
    #### GSA
    # GSE21340.22 -->Warning messages:
    # 1: In init.fit$sd < s0 :
    #   longer object length is not a multiple of shorter object length
    # 2: In sd + s0 :
    #   longer object length is not a multiple of shorter object length
    source('src/EAmethods/gsa.R')
    tmp_enrichment_df = 
      tryCatch(
        run_gsa(newSe, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = names(kegg_list)
      tmp_enrichment_df = tmp_enrichment_df %>%
        add_column(pathway_id = p,
                   GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,ES,Pvalue)
    } else {return()}
    
  } else if (method_name %in% c('gsva','replaid.gsva')){
    ##########################################
    #### GSVA
    source('src/EAmethods/gsva.R')
    use.plaid <- (method_name == 'replaid.gsva')
    tmp_enrichment_df = 
      tryCatch(
        run_gsva(newSe, kegg_list, matG, x[['seed_n']],
          use.plaid = use.plaid ),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = P.Value) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,t,Pvalue)
    } else {return()}
  } else if (grepl("^plaid",method_name)){
    ##########################################
    #### GSVA
    source('src/EAmethods/plaid.R')
    pars <- strsplit(method_name,split="[.]")[[1]]
    method="ttest"
    if(length(pars)>=2) method <- pars[2]
    tmp_enrichment_df = 
      tryCatch(
        run_plaid(newSe, matG, method=method),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue = pval) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,Pvalue)
    } else {return()}
    
  } else if (method_name == 'padog'){
    ##########################################
    #### PADOG
    source('src/EAmethods/padog.R')
    tmp_enrichment_df = 
      tryCatch(
        run_padog(newSe, kegg_list, 120, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      tmp_enrichment_df = tmp_enrichment_df %>%
        magrittr::set_rownames(NULL) %>%
        dplyr::rename(Pvalue = Ppadog,
                      pathway_id = ID) %>%
        add_column(GEO = paste0(x[['GEO']],'.',iter)) %>%
        mutate(meanAbsT0 = round(meanAbsT0,2),
               padog0 = round(padog0,4)) %>%
        select(GEO,pathway_id,meanAbsT0,padog0,PmeanAbsT,Pvalue)
    } else {return()}
    
  } else if (method_name == 'roast'){
    ##########################################
    #### ROAST
    source('src/EAmethods/roast.R')
    tmp_enrichment_df = 
      tryCatch(
        run_roast(newSe, kegg_list, x[['seed_n']]),
        error=function(cond) {
          return(NULL)})
    if (!is.null(tmp_enrichment_df)){
      p = rownames(tmp_enrichment_df)
      tmp_enrichment_df = tmp_enrichment_df %>%
        dplyr::rename(Pvalue=PValue) %>%
        magrittr::set_rownames(NULL) %>%
        add_column(pathway_id = p,
                   GEO = paste0(x[['GEO']],'.',iter)) %>%
        select(GEO,pathway_id,PropDown,PropUp,Direction,Pvalue)
    } else {return()}
  }
  
  ### FOR ALL METHODS
  return(tmp_enrichment_df)
}

enrichment_df = data.frame()
print(paste0('Datasets #',nrow(dataset_df)))
count = 1
tot = nrow(dataset_df)
d = dataset_df$GEO[1]
message("starting computation ",method_name)
for (d in dataset_df$GEO){
  message(paste0(count,'/',tot,' ', d)); count = count + 1 
  
  #################################################
  ## EXECUTE METHOD
  if(0) {
    pb = utils::txtProgressBar(min=0, max=30, style = 3)
    progress = function(n) utils::setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    cl = parallel::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)  
    boot = foreach::foreach(i = seq(30), .options.snow = opts)
    permuted_enrichment_df = foreach::`%dopar%`(boot, run_enrichment(d,i)) %>% 
      data.table::rbindlist()  
    parallel::stopCluster(cl)
    enrichment_df = rbind(enrichment_df, permuted_enrichment_df)
  } else {
    df <- run_enrichment(d, 1, method_name)
    enrichment_df = rbind(enrichment_df, df)
  }
}
enrichment_df %>%
  write_tsv(file = paste0(outdir,method_name,'.tsv.gz'))
message("finished!")
