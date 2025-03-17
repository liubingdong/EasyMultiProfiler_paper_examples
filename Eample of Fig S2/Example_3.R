library(EasyMultiProfiler)
library(magrittr)

# This example is to show the cached and history function.
# More detailed information are available on http://easymultiprofiler.xielab.net/

# load the meta data
data(MAE)

# 1.Figure S2A for the history tracing

k1 <- MAE %>%
  EMP_assay_extract('taxonomy') %>%
  EMP_collapse(estimate_group = 'Genus',collapse_by = 'row') %>%
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) %>%
  EMP_filter(feature_condition = pvalue<0.05)

k2 <- MAE %>%
  EMP_collapse(experiment = 'untarget_metabol',na_string=c('NA','null','','-'),
               estimate_group = 'MS2kegg',method = 'sum',collapse_by = 'row') %>%
  EMP_diff_analysis(method='DESeq2', .formula = ~Group) %>%
  EMP_filter(feature_condition = pvalue<0.05 & abs(fold_change) > 1.5)

#For two experinemnts 
p1 <- (k1 + k2) %>% EMP_cor_analysis(method = 'spearman') %>%
  EMP_heatmap_plot() ## Visualization

p1 %>% EMP_history()

# 2.Figure S2B for Caching Mechanism

## First run
system.time({
  MAE |>
    EMP_assay_extract('geno_ec') |>
    EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85)
})

## Second run
system.time({
  MAE |>
    EMP_assay_extract('geno_ec') |>
    EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85)
})
 
## Close the cached
system.time({
  MAE |>
    EMP_assay_extract('geno_ec') |>
    EMP_WGCNA_cluster_analysis(RsquaredCut = 0.85,use_cached = FALSE)
})
