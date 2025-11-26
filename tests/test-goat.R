install.packages("goat")
library(goat)

#source("~/Playground/playbase/dev/include.R",chdir=TRUE)
source("../R/ultragsea.R")
source("../R/gsetcor.R")
source("../R/gmt-utils.R")
source("../R/utils.R")

# TODO: change the output directory to some folder on your computer
# e.g. "~/data/goat" on unix systems, "C:/data/goat" on Windows,
# or set to getwd() to write output to the current working directory
output_dir = getwd()

# download an example gene list
datasets = download_goat_manuscript_data(output_dir)
genelist = datasets$`Wingo 2020:mass-spec:PMID32424284`

# download GO gene sets
genesets_asis = download_genesets_goatrepo(output_dir)

# filter gene sets for sufficient overlap with the gene list
genesets_filtered = filter_genesets(genesets_asis, genelist)

## genesets: tibble with genesets, must contain columns 'source',
##           'source_version', 'id', 'name', 'genes', 'ngenes',
##           'ngenes_signif'
head(genesets_filtered)
genesets_filtered$genes[[1]]

## genelist: tibble with genes, must contain column 'gene' and 'test'.
##           gene = character column, which are matched against list
##           column 'genes' in genesets tibble. test = boolean column (you
##           can set all to FALSE if not performing Fisher-exact or
##           hypergeometric test downstream)

head(genelist$gene)
head(genelist$test)
head(genelist$effectsize)


# apply GOAT; score_type indicates the input gene list should be ranked by respective
# "effectsize". This is the recommended setting. Alternatively, this can be set to "pvalue"
# as well but will yield fewer genesets (c.f. GOAT manuscript, real-world benchmarks)
result = test_genesets(genesets_filtered, genelist, method = "goat",
  score_type = "effectsize", padj_method = "bonferroni", padj_cutoff = 0.05)

# print first 10 rows of the result table
print(result |> select(source, name, ngenes, pvalue_adjust) |> utils::head(n=10))

# store results. This function will create 2 files:
# - an output table with gene set p-values
# - a log file that includes Methods text adapted to your settings
save_genesets(result, genelist, filename = paste0(output_dir, "/goat.xlsx"))


G <- Matrix::t(playdata::GSETxGENE)
fc <- Matrix::rowMeans(G!=0)
fc <- head(fc,8000)
gg <- intersect(names(fc),rownames(G))
fc <- fc[gg]
G <- G[gg,]
dim(G)
gmt <- mat2gmt(G)
gmt <- gmt[sapply(gmt,length) >= 10]

goat_test <- function(fc, gmt) {
  genesets  <- data.table::data.table(
    source = "GSET", source_version = "local",
    id = paste0("gset",1:length(gmt)),
    name = names(gmt),
    genes = gmt,
    ngenes = sapply(gmt,length),
    ngenes_signif = 0L
  )
  genesets <- as_tibble(genesets)
  
  genelist <- tibble(
    gene = names(fc),
    test = FALSE,
    signif = FALSE,  
    effectsize = fc
  )
  
  result <- test_genesets(
    genesets, genelist,
    method = "goat", score_type = "effectsize",
    padj_method = "bonferroni", padj_cutoff = 0.05)

  df <- data.frame(
    pathway = result$name,
    pval = result$pvalue,
    padj = result$pvalue_adjust,
    score = result$zscore,  
    size = result$ngenes
  )
  df
}

