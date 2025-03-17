library(EasyMultiProfiler)

# This data in the example is from  the microbiome-metabolome dataset from the Borenstein Lab
# https://github.com/borenstein-lab/microbiome-metabolome-curated-data

## We have prepared the data in the standard MultiAssayExperment format
Bowel <- readRDS('Bowel.rds')

## Perform the microbial analysis
## Due to the big-scale data, this process may cost about 3~6 minutes
Bowel |>
  EMP_assay_extract('taxonomy') |>
  EMP_adjust_abundance(method = 'combat',
                       .factor_unwanted = 'Dataset',
                       .factor_of_interest = 'Study.Group') |>  ### Remove studies batch effect
  EMP_decostand(method = 'relative') |>
  EMP_collapse(estimate_group = 'Genus',
               collapse_by = 'row') |>
  EMP_alpha_analysis() |> 
  EMP_boxplot(estimate_group = 'Study.Group',plot='violin',
              ref.group = 'Control') |>
  EMP_dimension_analysis(method = 'pcoa',
                         distance = 'bray') |>
  EMP_scatterplot(show = 'p12html',
                  force_adonis = TRUE,
                  ref.group = 'Control') |> ### The Adonis analysis cost much time
  EMP_identify_assay(method = 'default') |>  ### Filter the noise taxa
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

## Perform the metabolitecs analysis
Bowel |>
  EMP_assay_extract('metabolite') |> 
  EMP_filter(Study.Group %in% c('Control','CRC')) |> ### Select the Control and CRC
  EMP_collapse(estimate_group = 'KEGG',collapse_by = 'row') |>
  EMP_adjust_abundance(method = 'combat',
                       .factor_unwanted = 'Dataset', 
                       .factor_of_interest = 'Study.Group') |> ### Remove studies batch effect
  EMP_decostand(method = 'relative') |>
  EMP_WGCNA_cluster_analysis(deepSplit = 4,
                             minModuleSize = 10,
                             mergeCutHeight = 0.1) |> ### WGCNA cluster
  EMP_save_var('CRC_metabolite') |>
  EMP_identify_assay(method = 'default') |>
  EMP_diff_analysis(method = 'wilcox.test',
                    estimate_group = 'Study.Group') |>  
  EMP_volcanol_plot(dot_size = 4) |>
  EMP_dimension_analysis(method = 'opls') |> ### OPLS model
  EMP_scatterplot(show = 'p12html',
                  force_adonis = TRUE) |>  
  EMP_filter(feature_condition = VIP >1) |>  ### Select the marker metabolite
  EMP_enrich_analysis(keyType = 'cpd',pvalueCutoff = 0.05,
                      KEGG_Type = 'KEGG') |> ### KEGG enrichment
  EMP_enrich_dotplot(show = 10) 

(CRC_metabolite + all_tax |> ### Combime microbiome and metabolite data
    EMP_filter(Study.Group %in% c('Control','CRC')) ) |>
  EMP_WGCNA_cor_analysis() |> ### WGCNA cor analysis
  EMP_heatmap_plot(label_size = 2) |>
  EMP_assay_extract('metabolite') |>
  EMP_filter(feature_condition = WGCNA_color %in% c('blue',
                                                    'brown',
                                                    'turquoise')) |> ### Select the hub
  EMP_enrich_analysis(keyType = 'cpd',
                      pvalueCutoff = 0.05,
                      KEGG_Type = 'MKEGG') |>
  EMP_enrich_dotplot()


