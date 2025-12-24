
methods = c(
#  'ultragsea.cor',
  'ultragsea.cor',
  'ultragsea.ztest', 
  'cameraPR','goat', 'fgsea','fisher', ##'gsva',
  'fastFET', 'plaid', 'replaid.gsva')

methods = sort(methods)
m="goat"
m="fastFET"
m="fisher"
m="cameraPR"
m="ultragsea.cor"
m="ultragsea.ztest"

model_name=m

for(m in methods) {
  if(m %in% c("fisher","fastFET")) {
    system(paste("Rscript src/EAcategory/TPbenchmark/OVA.R -m",m))
    system(paste("Rscript src/EAcategory/FPbenchmark/OVA.R -m",m))
  } else {
    system(paste("Rscript src/EAcategory/TPbenchmark/PGA.R -m",m))
    system(paste("Rscript src/EAcategory/FPbenchmark/PGA.R -m",m))
  }
}

# This script essentially consolidates all enrichment analysis
# predictions into standardized formats for subsequent performance
# evaluation.
source("src/PerformanceEvaluation/gather_results.R")
source("src/PerformanceEvaluation/extract_tpr&fpr.R")
source("src/PerformanceEvaluation/extract_ranking.R")

## make plots
source("notebooks/runBenchmarkFiguresAndTables.R")

