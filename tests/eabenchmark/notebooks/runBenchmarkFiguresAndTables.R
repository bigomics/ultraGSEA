rm(list = setdiff(ls(),"methods"))

library(tidyverse)
library(knitr)
library(scales)
library(reshape2)
library(ggridges)

reorder_within=function(x, by, within, fun=median, sep='___', ...) {
  new_x=paste(x, within, sep=sep)
  stats::reorder(new_x, by, FUN=fun)
}

scale_x_reordered=function(..., sep='___') {
  reg=paste0(sep, '.+$')
  ggplot2::scale_x_discrete(labels=function(x) gsub(reg, '', x), ...)
}

symlog_trans <- function(base = 10, thr = 1, scale = 1){
  ##This is an R function that creates a symmetric logarithmic
  ## transformation (symlog_trans) for use with plotting data in a
  ## graph.
  trans <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             (thr + scale * suppressWarnings(log(sign(x) * x / thr, base))))

  inv <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             base^((sign(x) * x - thr) / scale) * thr)

  breaks <- function(x){
    sgn <- sign(x[which.max(abs(x))])
    if(all(abs(x) < thr))
      pretty_breaks()(x)
    else if(prod(x) >= 0){
      if(min(abs(x)) < thr)
        sgn * unique(c(pretty_breaks()(c(min(abs(x)), thr)),
                       log_breaks(base)(c(max(abs(x)), thr))))
      else
        sgn * log_breaks(base)(sgn * x)
    } else {
      if(min(abs(x)) < thr)
        unique(c(sgn * log_breaks()(c(max(abs(x)), thr)),
                 pretty_breaks()(c(sgn * thr, x[which.min(abs(x))]))))
      else
        unique(c(-log_breaks(base)(c(thr, -x[1])),
                 pretty_breaks()(c(-thr, thr)),
                 log_breaks(base)(c(thr, x[2]))))
    }
  }
  trans_new(paste("symlog", thr, base, scale, sep = "-"), trans, inv, breaks)
}


## To test your own method:
##  1. Add your method name to methods variable
##  2. Add your method to the assign_category function
##  3. Add your method name to the convert_name function
##  4. Edit plot function in cell SkewenessHistogram to accept 1 more row (facet_wrap with nrow>5)

networks=c('funcoup','string')
ora = c('ease','fisher')
anubix = paste('anubix',networks,sep = '-')
neat = paste('neat',networks,sep = '-')
binox = paste('binox',networks,sep = '-')
netpea = paste('netpea',networks,sep = '-')
camera = c('camera','camera*')
#methods = c(anubix,binox,neat,netpea,camera,ora,'roast','fgsea','gsa','gsea','padog','gsva','cepaORA','spia')
#methods = c('cameraPR','goat','ultragsea')

assign_category = function(m){
  if(m == 'fisher') return('Overlap-based')
  if(m == 'fastFET') return('Overlap-based')
  if(m == 'fgsea') return('Preranked (PR)')
  if(m == 'goat') return('Preranked (PR)')
  if(m == 'cameraPR') return('Preranked (PR)')
  if(m == 'gsva') return('Single-sample method (SS)')
  if(grepl("ultragsea",m)) return('Preranked (PR)')  
  if(grepl('plaid',m)) return('Single-sample method (SS)')
}

convert_name = function(m){
  if(m == 'fisher') return('Fisher')
  if(m == 'fgsea') return('fGSEA')
  if(m == 'goat') return('GOAT')
  if(m == 'cameraPR') return('cameraPR')
  if(m == 'gsva') return('GSVA')
  if(m == 'replaid.gsva') return('replaidGSVA')  
  #if(m == 'plaid') return('PLAID')
  if(m == 'ultragsea.cor') return('ultraGSEA.cor')
  if(m == 'ultragsea.ztest') return('ultraGSEA.ztest')
  return(m)
}

## tool_shape=c(
##   'fGSEA'=0,
##   'GOAT'=21,
##   'cameraPR'=22,
##   'GSVA'=23,
##   'PLAID'=24,
##   'replaidGSVA'=25,
##   'Fisher'=21,
##   'ultraGSEA.cor'=1,
##   'ultraGSEA.ztest'=2
## )
tool_shape <- head( rep(c(0:6,21:25,7:20),99),length(methods))
names(tool_shape) <- sapply(methods,convert_name)
tool_shape = tool_shape[!duplicated(names(tool_shape))]
tool_shape = tool_shape[order(names(tool_shape))]

assignCategory <- function(m) {
  m2cat <- Vectorize(assign_category)(unique(m))
  factor( m2cat[m],
    levels=c('Overlap-based',
      'Preranked (PR)',
      'Single-sample method (SS)',
      'FunCoup-based method',
      'STRING-based method',
      'Pathway Topology-based')
  )
} 
convertName <- function(m) {
  new.name <- Vectorize(convert_name)(unique(m))
  factor( new.name[m], levels=names(tool_shape))
} 


## color palette
# https://coolors.co/bee9e8-62b6cb-1b4965-cae9ff-5fa8d3
category_col= c(
   'Overlap-based'= '#588157',
   'STRING-based method'= '#ff8fab',
   'FunCoup-based method'= '#014f86',
   'Pathway Topology-based'='#f9dc5c',
   'Preranked (PR)'='#6095bb',
   'Single-sample method (SS)'='#e63946'
)

# Define a plotting theme for use throughout the document
plain_gg_theme=theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.border = element_rect(linetype = "solid", colour = "black", size=0),
          axis.line = element_line(colour = 'black', size = 1),
          axis.ticks = element_line(colour = "black", size = 1),
          plot.title=element_text(size=13, face='bold'),
          axis.title=element_text(size=13, face='bold'),
          strip.text=element_text(size=13, face='bold'),
          strip.background = element_rect(colour="white", fill="white"),
          legend.justification=c('center', 'top'),
          legend.box.just='right',
          legend.margin=margin(6, 6, 6, 6),
          legend.title= element_text(size=10, face='bold'),
          legend.text= element_text(size=10))


## Benchmark datasets

## We retrieved and curated 82 high-quality datasets from
## [Gemma](https://gemma.msl.ubc.ca), focusing on case-control
## observational studies. These datasets, sourced from GEO, underwent
## rigorous curation, ensuring reliable gene expression data. Our
## criteria included a minimum of three samples per condition, no drug
## treatment, batch effect correction, and at least 15 differentially
## expressed genes (FDR < 0.2). We processed DNA-microarray data through
## quantile normalization and log transformation, while RNA-Seq data were
## obtained as log2-transformed counts per million reads. We mapped
## probes to Gene IDs, and for genes with multiple probes, we averaged
## the values. We performed quality checks and analyses, detailed in
## Supplementary Materials.

#################################################
## EXPRESSION DATA
dataset_df = 
  read_tsv(file='data/TPbenchmark/cc_deg-alpha10OR20-Min15Max500-beta0-DEG.tsv', 
           col_types = 'fcfffiddiciccccc') %>%
  arrange(TargetPathway)

dataset_df %>%
  mutate(NrOfSamples = paste0(NrOfSamples, ' (', CC,')'),
         NrOfDEG = paste0(NrOfDEG, ' (', UpDown,')')) %>%
   dplyr::select(GEO,TargetPathway,NrOfSamples,GenomeCoverage,NrOfDEG,BatchEffect) %>%
  mutate(BatchEffect = str_replace(BatchEffect, "BATCH_CORRECTED_SUCCESS|Data has been batch-corrected", "Corrected")) %>%
  mutate(BatchEffect = str_replace(BatchEffect, "No batch effect was detected", "Not detected")) %>%
   dplyr::rename(`Batch effect` = BatchEffect, 
         `Disease/Target pathway` = TargetPathway,
         `Nr. of samples`=NrOfSamples,
         `Nr. of DEGs`=NrOfDEG,
         `Genome coverage` = GenomeCoverage) %>% 
  write_tsv('notebooks/tables/GeneExpressionMetadata.tsv') %>% 
  knitr::kable(format = 'simple') 


## Disease Pathway Network

## We introduce the Disease Pathway Network to counteract the shortcoming
## of the "single target pathway" approach and improve the sensitivity
## assessment in an unbiased way. HumanNet-XC, a comprehensive functional
## network of human genes, was used to analyze inter-pathway
## connectivity. Inter-pathway connectivity (IPC) was quantified as the
## sum of direct links and shared neighbors between genes in two
## pathways, tested for significance using subsampling. Inter-pathway
## overlap was assessed using the Jaccard index. Pathway pairs with BH
## FDR-corrected p-values < 0.05 in both tests were retained.

kegg_df = read_tsv(file = 'data/KEGG/keggrest_hsa.tsv', col_types = 'ffcccff')
kegg_df = kegg_df %>%
  filter(pathway_id %in%
           as.character(
             kegg_df %>%
                dplyr::select(pathway_id,entrezgene_id) %>%
               unique() %>%
               group_by(pathway_id) %>%
               summarise(n=n()) %>%
               filter(n>14) %>%
               droplevels() %>%
               pull(pathway_id))) %>%
  droplevels() %>%
  dplyr::rename(GeneID=entrezgene_id)

kegg_list = tapply(kegg_df$GeneID, kegg_df$pathway_id, function(x) as.character(x))
kegg_list = lapply(kegg_list, function(x) unique(x))
kegg_list = kegg_list[sapply(kegg_list,length) >= 10]

kegg_df = kegg_df %>% 
  dplyr::select(pathway_id,pathway_name,class,subclass) %>% 
  unique()

kegg_class_df = kegg_df %>%
  dplyr::select(pathway_id,pathway_name,class,subclass) %>% 
  dplyr::rename(PathwayClass = 'class',
         PathwaySubclass = 'subclass') %>%
  unique()


top_cutoff = c(10,20,40,'ALL')
tmp_pathway_network_df = lapply(top_cutoff, function(top){
    tmp_df = read_tsv(file = paste0('results/DiseasePathwayNetwork/DiseasePathwayNetwork-FDR<5e-2_top',top,'GOaverageBMA.tsv.gz'), col_types = cols())
  tmp_df$TOP = top
  return(tmp_df)
}) %>% data.table::rbindlist()

pathway_network_df = tmp_pathway_network_df %>% 
  group_by(A_name, TOP) %>%
  summarise(n = n()+1) %>% # +1 is for Target Pathway
  ungroup() %>%
  pivot_wider(names_from=TOP, values_from=n) %>% 
  dplyr::rename(pathway_name = A_name) %>% 
  left_join(kegg_class_df, by='pathway_name') %>%
  replace(is.na(.), 1) %>% 
  dplyr::select(pathway_name,pathway_id,PathwaySubclass,`10`,`20`,`40`,'ALL') %>%
  dplyr::rename('Pathway name' = pathway_name,
         'Pathway ID' = pathway_id,
         'Pathway subclass (Human diseases)' = PathwaySubclass)

write_tsv(pathway_network_df,'notebooks/tables/DiseasePathwayNetwork.tsv')

pathway_network_df %>%
  dplyr::filter(`Pathway name` %in% unique(dataset_df$TargetPathway)) %>%
  write_tsv('notebooks/tables/TargetDiseasePathwayNetwork.tsv') %>% 
  knitr::kable(format = 'simple') 


##========================================================
## Nr. of enrichments per dataset
##========================================================

## Min, max and median average number of tested pathways in the positive benchmark.

#positive_benchmark_count = read_tsv(file = 'notebooks/tables/TPbenchmark_counts.tsv', col_types = 'cif')
positive_benchmark_count = read_tsv(file = 'results/stats/TPbenchmark_counts.tsv', col_types = 'cif')
order_m = positive_benchmark_count %>% filter(type=='median') %>% arrange(-val) %>% pull(Method)
order_m = as.character(sapply(order_m, function(x) convert_name(x)))
positive_benchmark_count$Method=sapply(positive_benchmark_count$Method, function(x) convert_name(x))
positive_benchmark_count$Method = factor(positive_benchmark_count$Method, levels = order_m)
table(positive_benchmark_count$Method)
head(positive_benchmark_count)

pdf("plots/NrEnrichment.pdf",w=10,h=5)
ggplot(positive_benchmark_count, aes(x=Method, y=val, fill = type)) +
  geom_bar(position="dodge", stat="identity",width = 0.7, color = "black")+
  xlab('')+
  ylab('Number of tested pathways') +
  scale_fill_manual(values = c('min'='#26547c', 'max'='#ef476f', 'median'='#ffd166')) +
  plain_gg_theme + 
  theme(axis.title.y=element_text(size=13,face='bold',colour="black"),
        axis.text.y=element_text(size=13,face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold', angle=45, hjust =1, colour="black"),
        legend.title = element_blank(),
        plot.margin=unit(c(1,1,0,1), 'cm'))  
dev.off()

##-------------------------------------------------------
## Analysis of biases in the benchmarked methods
##-------------------------------------------------------

## We assessed the EA methods' performance on randomized data. In an
## ideal scenario, the method should produce a uniform distribution of
## p-values ranging from 0 to 1 for pathways when applied to
## randomized data. Ideally, 5% of these p-values would fall below the
## cutoff of 0.05.
balanced = TRUE
if(balanced) {
  negative_benchmark=read_tsv('results/stats/FPbenchmark_balanced.tsv.gz') %>%
    replace(is.na(.),1) %>%
    filter(Method %in% methods) %>% 
    droplevels()
} else {
  negative_benchmark=read_tsv('results/stats/FPbenchmark.tsv.gz') %>%
    replace(is.na(.),1) %>%
    filter(Method %in% methods) %>% 
    droplevels()
}

## check the fraction of pvalue below 5% per pathway
df=negative_benchmark  %>%
  group_by(Method,pathway_id) %>%
  summarise(FPR=sum(Pvalue<0.05)/n()) %>%
  ungroup()

#df$Category = sapply(df$Method, function(x) assign_category(x))
df$Category = assignCategory(df$Method)
#df$Method = sapply(df$Method, function(x) convert_name(x))
df$Method = convertName(df$Method)

# df %>%
#   group_by(Method) %>%
#   summarise(round(median(FPR),3))

pdf("plots/PvalueUnderTrueH0.pdf",w=6,h=4)
ggplot(df, aes(x=reorder_within(Method, by=FPR, within=Method, fun=median),y=FPR, fill=Category)) +
  geom_boxplot(outlier.size=0.5) +
  scale_x_reordered() +
  xlab('') + 
  scale_y_continuous("False Positive Rate",
                   breaks = seq(0,1,0.1)) +
  scale_fill_manual(values=category_col) +
  geom_hline(aes(yintercept=0.05), linetype='dashed',color='black', alpha =0.5) +#
  plain_gg_theme +
  theme(axis.title.y=element_text(size=13,face='bold',colour="black"),
        axis.text.y=element_text(size=13,face='bold',colour="black"),
        axis.text.x=element_text(size=12, face='bold', angle=45, hjust =1, colour="black"),
        legend.position='none',
        plot.margin=unit(c(1,1,0,1), 'cm'))
dev.off()


## Under the null hypothesis, EA methods often yield p-values that
## display a bias either toward 0 or 1 or exhibit a bimodal distribution
## skewed towards both extremes. This bias can significantly influence
## the significance of the analysis. Therefore, we studied the p-value
## distributions for each method to determine if they were right-skewed
## (biased toward 0) or left-skewed (biased toward 1). A right-skewed
## distribution (p-values biased toward 0) has the potential to produce
## false positives by identifying pathways as affected when they are
## not. Conversely, a left-skewed distribution (p-values biased toward 1)
## may lead to false negatives by indicating pathways as non-significant
## when they are indeed impacted.


#### Visualize SKEWNESS OF PATHWAY NULL HYPO
df=negative_benchmark %>%
  select(Pvalue,Method) %>%
  na.omit()

df$Category = assignCategory(df$Method)
df$Method = convertName(df$Method)
## df$Category=sapply(df$Method, function(x) assign_category(x))
## df$Method=sapply(df$Method, function(x) convert_name(x))
## df$Category=
##   factor(df$Category, 
##          levels=c('Overlap-based',
##                   'Preranked (PR)',
##                   'Single-sample method (SS)',
##                   'FunCoup-based method',
##                   'STRING-based method',
##                   'Pathway Topology-based'))
## df$Method=factor(df$Method, levels=names(tool_shape))


df$Pvalue[df$Pvalue==1] <- NA

pdf("plots/SkewenessHistogram.pdf", w=12, h=8)
ggplot(df, aes(x = Pvalue, fill = Category)) +
    geom_density(stat = "density",kernel='gaussian') +
    facet_wrap(~Method, scales = 'free_y', nrow = 3, ncol = 4) + 
    scale_fill_manual(values=category_col) +
    ylab('Density') +
    xlab('P-value') +
    plain_gg_theme + 
    theme(axis.title.x=element_text(size=15,face='bold',colour="black"),
          axis.title.y=element_text(size=15,face='bold',colour="black"),
          axis.text.y=element_text(size=13,face='bold',colour="black"),
          axis.text.x=element_text(size=13, face='bold',colour="black"),
          strip.text=element_text(size=15, face='bold'),
          legend.position='none',
          plot.margin=unit(c(1,1,1,1), 'cm'))
dev.off()


##======================================================================
## Stats at different number of TPs in the Disease Pathway Network
##======================================================================

## In this study, we used independent positive and negative
## benchmarks. The positive benchmark includes true positives (TP) and
## false negatives (FN), representing pathways correctly identified as
## significant (p-value < 0.05) and non-significant (p-value ≥ 0.05),
## respectively. Similarly, the negative benchmark includes true
## negatives (TN) and false positives (FP), indicating pathways correctly
## identified as non-significant or significant, respectively. We created
## the negative benchmark by resampling gene labels on the datasets from
## the genome, ensuring a consistent number of differentially expressed
## genes (DEGs) for accurate false positive rate (FPR) estimation across
## tests. To address the imbalance between positive and negative
## pathways, we focused on target-related pathways in the negative
## benchmark to calculate TN and FP.

## Using TP, TN, FP, and FN definitions, we derived true positive rate
## (TPR, or sensitivity) and true negative rate (TNR, or specificity, or
## 1-FPR). We computed the geometric mean of TPR and TNR (G-mean) as a
## comprehensive performance summary. Additionally, we assessed the
## median relative rank of TPs among the top predictions, considering
## ties by averaging the ranks.

top_cutoff=c(10,20,40,'ALL')

stats_df=read_tsv('results/stats/stats_balanced.tsv.gz') %>%
  filter(Method %in% methods) %>%
  replace(is.na(.), 0)
# 
rank_df=read_tsv('results/stats/Ranking(AVG).tsv.gz') %>%
  filter(Method %in% methods) %>%
  replace(is.na(.), 1)


medianRank_df=rank_df %>% 
              group_by(Method,TOP) %>%
              summarise(mRK=median(rank_all_mean, na.rm=T)) %>%
              ungroup()

tmp_total_stats_df=stats_df %>%
  left_join(medianRank_df, by=c('Method','TOP')) %>%
  mutate(mRK = round(mRK,2),
         GM = round(GM,2))
tmp_total_stats_df$Category=sapply(tmp_total_stats_df$Method, function(x) assign_category(x))
tmp_total_stats_df$Method=sapply(tmp_total_stats_df$Method, function(x) convert_name(x))

tmp_total_stats_df$Category=
  factor(tmp_total_stats_df$Category, 
         levels=c('Overlap-based',
                    'Preranked (PR)',
                    'Single-sample method (SS)',
                    'FunCoup-based method',
                    'STRING-based method',
                    'Pathway Topology-based'))

tmp_total_stats_df$TOP=ifelse(tmp_total_stats_df$TOP=='ALL', 'All links', paste0('Top ',tmp_total_stats_df$TOP, ' links'))
tmp_total_stats_df$TOP = factor(tmp_total_stats_df$TOP, levels=c('Top 10 links','Top 20 links','Top 40 links','All links'))
tmp_total_stats_df$Method=as.factor(tmp_total_stats_df$Method)
tmp_total_stats_df$Category=as.factor(tmp_total_stats_df$Category)

category_col_df = 
  data.frame(Category = names(category_col),
           colour = as.character(category_col))

total_stats_df = tmp_total_stats_df %>%
  left_join(category_col_df,by = 'Category') %>%
  mutate(size = 4,stroke = 2)

tool_colour = total_stats_df$colour
names(tool_colour) = total_stats_df$Method
tool_frame_colour = total_stats_df$colour
names(tool_frame_colour) = total_stats_df$Method
tool_size = total_stats_df$size
names(tool_size) = total_stats_df$Method

tool_colour = tool_colour[unique(names(tool_colour))]
tool_size = tool_size[unique(names(tool_size))]
total_stats_df$Method = factor(total_stats_df$Method, levels=names(tool_shape))


head(stats_df)
stats_dd <- data.frame(stats_df)
stats_dd <- apply(stats_dd[,-c(1,2)], 2, function(x) tapply(x, stats_dd$Method, mean))
stats_dd <- data.frame(stats_dd)
head(stats_dd)
stats_dd$Method <- sapply(rownames(stats_dd), function(x) convert_name(x))
stats_dd$Category <- sapply(rownames(stats_dd), function(x) assign_category(x))

pdf("plots/TPRvsFPR.pdf", height=6, width=5.5)
ggplot(stats_dd, aes(x=TPR, y=1-FPR, shape=Method, fill=Method, color=Method,
  size=Method, stroke=2 )) + 
  geom_point() +
  scale_shape_manual(values = tool_shape) +
  scale_fill_manual(values = tool_colour) +
  scale_color_manual(values=tool_frame_colour) +
  scale_size_manual(values=tool_size) +
  guides(color = guide_legend(override.aes = list(stroke = 2, size = 2.5))) + 
  xlab('TPR') + 
  ylab('1 - FPR') +
  plain_gg_theme +
  theme(axis.text.y=element_text(size=13,face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold',colour="black"),
        legend.position = "top", 
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=9, face='bold',colour="black"))
dev.off()


##======================================================================
## Runtime
##======================================================================

## We conducted scalability tests for each method in the benchmark using
## KEGG as input on 82 datasets. Some methods support parallelization and
## can handle multiple datasets simultaneously, reducing elapsed time in
## a battery testing setup. The analysis was conducted on macOS Monterey
## (v.12.5.1) with an Apple M1 processor (16GB RAM), except for BinoX,
## which ran on Ubuntu (v.18.04.6) with an Intel Core i7-2600 3.40GHz
## processor (16GB RAM). GSEA is presented as elapsed time in the
## results.

### Runtime
runtime_df=lapply(methods,function(m){
  df=read_tsv(paste0('results/runtime/',m,'.tsv.gz'), col_types='fddd')
  df$Method=m  
  if(m=='gsea'){
   df = df %>%
      dplyr::select(GEO,elapsed,Method) %>%
      dplyr::rename(user=elapsed)
  }
  return(df %>% 
           dplyr::select(GEO,user,Method))
}) %>% data.table::rbindlist()

runtime_df$Category=sapply(runtime_df$Method, function(x) assign_category(x))
runtime_df$Method=sapply(runtime_df$Method, function(x) convert_name(x))

order_runtime_df = runtime_df %>%
  group_by(Method) %>%
  summarise(muser = median(user,na.rm = T)) %>% 
  arrange(muser)

runtime_df$Method = factor(runtime_df$Method, levels = order_runtime_df$Method)
sort(tapply(runtime_df$user, runtime_df$Method, mean, na.rm=TRUE))
#runtime_df$user <- runtime_df$user + 1e-8
#runtime_df$user[runtime_df$user < 1e-3] <- NA

pdf("plots/RuntimeDataset.pdf", height=4, width=7)
ggplot(runtime_df, aes(x=user, y=Method, fill=Category)) +
  geom_density_ridges() +
  scale_x_continuous("CPU time (s)",
    #transform = symlog_trans(),
    transform = "log10"
    #breaks = c(0,0.2,0.5,1,10,100,1000)
    ) +
  scale_fill_manual(values=category_col) +
  ylab('') +
  # ggtitle('OBS: GSEA is shown as elapsed time') + 
  plain_gg_theme + 
  theme(axis.title.y=element_text(size=13,face='bold',colour="black", vjust=5),
        axis.text.y=element_text(size=12,face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold',colour="black"),
        legend.position='none',
        plot.margin=unit(c(1,1,1,1), 'cm'))
dev.off()


avg.runtime <- tapply(runtime_df$user, runtime_df$Method, mean, na.rm=TRUE)
stats_dd$avg.runtime <- avg.runtime[stats_dd$Method]
stats_dd$speedup <- mean(stats_dd$avg.runtime,na.rm=TRUE) / stats_dd$avg.runtime

pdf("plots/TPRvsSpeed.pdf", height=6, width=5.5)
ggplot(stats_dd,
  aes(x=speedup, y=TPR, shape=Method, fill=Method, color=Method,
  size=Method, stroke=2 )) + 
  geom_point() +
  scale_x_log10() +
  scale_shape_manual(values = tool_shape) +
  scale_fill_manual(values = tool_colour) +
  scale_color_manual(values=tool_frame_colour) +
  scale_size_manual(values=tool_size) +
  guides(color = guide_legend(override.aes = list(stroke = 2, size = 2.5))) + 
  xlab('relative speedup') + 
  ylab('TPR') +
  plain_gg_theme +
  theme(axis.text.y=element_text(size=13,face='bold',colour="black"),
        axis.text.x=element_text(size=13, face='bold',colour="black"),
        legend.position = "top", 
        legend.box = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=9, face='bold',colour="black"))
dev.off()
