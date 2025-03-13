library(EasyMultiProfiler)

Bowel <- readRDS('Bowel.rds')
system.time({
  

Bowel |>
  EMP_assay_extract('taxonomy') |>
  EMP_adjust_abundance(method = 'combat',
                       .factor_unwanted = 'Dataset',
                       .factor_of_interest = 'Study.Group') |>
  EMP_decostand(method = 'relative') |>
  EMP_collapse(estimate_group = 'Genus',
               collapse_by = 'row') |>
  EMP_alpha_analysis() |>
  EMP_boxplot(estimate_group = 'Study.Group',
              ref.group = 'Control') |>
  EMP_dimension_analysis(method = 'pcoa',
                         distance = 'bray') |>
  EMP_scatterplot(show = 'p12html',
                  force_adonis = F,
                  ref.group = 'Control') |>
  EMP_identify_assay(method = 'default') |>
  EMP_marker_analysis(method = 'boruta') |>
  EMP_filter(feature_condition = Boruta_decision!='Rejected') |>
  EMP_save_var('all_tax') |>
  EMP_collapse(estimate_group = 'Study.Group',
               method = 'mean',
               collapse_by = 'col') |>
  EMP_heatmap_plot(scale = 'standardize',
                   clust_row = TRUE,
                   clust_col = TRUE,
                   palette = 'Spectral') 
})

Bowel |>
  EMP_assay_extract('metabolite') |> 
  EMP_filter(Study.Group %in% c('Control','CRC')) |>
  EMP_collapse(estimate_group = 'KEGG',collapse_by = 'row') |>
  EMP_adjust_abundance(method = 'combat',
                       .factor_unwanted = 'Dataset', 
                       .factor_of_interest = 'Study.Group') |>
  EMP_decostand(method = 'relative') |>
  EMP_WGCNA_cluster_analysis(deepSplit = 4,
                             minModuleSize = 10,
                             mergeCutHeight = 0.1) |>
  EMP_save_var('CRC_metabolite') |>
  EMP_identify_assay(method = 'default') |>
  EMP_diff_analysis(method = 'wilcox.test',
                    estimate_group = 'Study.Group') |>
  EMP_volcanol_plot(dot_size = 4) |>
  EMP_dimension_analysis(method = 'opls') |>
  EMP_scatterplot(show = 'p12html',
                  force_adonis = TRUE) |>
  EMP_filter(feature_condition = VIP >1) |>
  EMP_enrich_analysis(keyType = 'cpd',pvalueCutoff = 0.05,
                      KEGG_Type = 'KEGG') |>
  EMP_enrich_dotplot(show = 10) 

(CRC_metabolite + all_tax |>
    EMP_filter(Study.Group %in% c('Control','CRC')) ) |>
  EMP_WGCNA_cor_analysis(method = 'pearson') |>
  EMP_heatmap_plot(label_size = 2) |>
  EMP_assay_extract('metabolite') |>
  EMP_filter(feature_condition = WGCNA_color %in% c('blue',
                                                    'brown',
                                                    'turquoise')) |>
  EMP_enrich_analysis(keyType = 'cpd',
                      pvalueCutoff = 0.05,
                      KEGG_Type = 'MKEGG') |>
  EMP_enrich_dotplot()


