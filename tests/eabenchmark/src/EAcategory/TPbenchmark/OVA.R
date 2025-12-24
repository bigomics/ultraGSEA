rm(list = ls())
setwd('~/eabenchmark/') 

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(parallel))

option_list = list(
  make_option(c("-c", "--ncores"), type="integer", default='1',
              help=" (by default: %default)", metavar="character"),
  make_option(c("-m", "--method"), type="character", default='ease',
              help=" (by default: %default)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
method_name = opt$method
ncores = opt$ncores

#################################################
## FUNTIONAL GENE SETS
message('Load KEGG')
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

matG <- ultragsea::gmt2mat(kegg_list)

#################################################
## EXPRESSION DATA
dataset_df = read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', col_types = 'fcfffiddiciccccc')
head(dataset_df)

outdir = paste0('results/TPbenchmark/')

#################################################
## EXECUTE METHOD

d=dataset_df$GEO[1]
d

run_enrichment = function(d){
  library(tidyverse)

  x = dataset_df[dataset_df$GEO==d,]

  se = readRDS(paste0('data/TPbenchmark/',x[['DataType']],'/',x[['GEO']],'.Rds'))
  bg <- rownames(se)
  genes <- as.character(unlist(stringr::str_split(x[['GeneID']],',')))

  kegg_list2 <- lapply(kegg_list, function(s) intersect(s,bg))
  genes <- intersect(genes, bg)
  
  if (method_name == 'ease'){
    source('src/EAmethods/ease.R')
    ptm = proc.time()
    tmp_enrichment_df = run_ease( genes, 
                                 PT=as.numeric(x[['GenomeCoverage']]),
                                 kegg_list2)
    tmp_enrichment_df$GEO = x[['GEO']]
  } else if (method_name == 'fisher'){
    source('src/EAmethods/fisher.R')
    ptm = proc.time()
    tmp_enrichment_df = run_fisher( genes, 
                                   PT=as.numeric(x[['GenomeCoverage']]),
                                   kegg_list2)
    tmp_enrichment_df$GEO = x[['GEO']]
  } else if (method_name == 'fastFET'){
    source('src/EAmethods/fisher.R')
    ptm = proc.time()
    gg <- which(rownames(matG) %in% bg)
    tmp_enrichment_df = run_fastFET(genes, matG[gg,], bg) %>%
      dplyr::rename(Pvalue = p.value) %>%
      magrittr::set_rownames(NULL) %>%
      add_column( pathway_id = colnames(matG),
        GEO = x[['GEO']] )
  } else if (method_name == 'random'){
    source('src/EAmethods/random.R')
    ptm = proc.time()
    tmp_enrichment_df = run_random(genes,
                                   kegg_list2)
    tmp_enrichment_df$GEO = x[['GEO']]
  }
  
  ### FOR ALL METHODS
  runtime = proc.time() - ptm
  tmp_enrichment_df$user=round(runtime[[1]],2)
  tmp_enrichment_df$system=round(runtime[[2]],2)
  tmp_enrichment_df$elapsed=round(runtime[[3]],2)
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

# output runtime
output_df %>%
  select(GEO,user,system,elapsed) %>%
  unique() %>%
  write_tsv(file = paste0('results/runtime/',method_name,'.tsv.gz'))
