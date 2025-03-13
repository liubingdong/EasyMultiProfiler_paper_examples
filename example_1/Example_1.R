library(EasyMultiProfiler)

meta_data <- read.csv('meta.csv',header = TRUE,row.names = 1,na.strings = "") 

MAE <- EMP_easy_import('GGMP7009_even10k.biom',file_format = 'biom',
                       coldata = meta_data,type = 'tax') 
MAE |> EMP_summary()

patient_physi <- c('BMI','FBG','HDL','LDL','TCHO','TG','UA')
patient_diet <- c('Beer','Soft_drinks','Animal_oil','Plant_oil',
                  'Salt','Sauce','Soy_sauce','Sugar','Fruit_juice',
                  'Fruits','Grains','H_alcohol','Livestock','L_alcohol',
                  'Rice_wine','Vegetables','Red_wine','Yellow_wine')
MAE |>
  EMP_assay_extract('experiment') |>
  EMP_coldata_extract(coldata_to_assay = patient_physi,action = 'add') |>
  EMP_save_var('meta_data') |>
  EMP_coldata_extract(coldata_to_assay = patient_diet,action = 'add') |>
  EMP_save_var('diet_data') |>
  EMP_cluster_analysis(h = 0.5) |>
  EMP_filter(cluster == 24) |>
  EMP_filter(!is.na(Metabolic_syndrome)) |>
  EMP_filter(Antibiotics == 'n' & Colitis == 'n') |>
  EMP_filter(Diarrhea =='n' & Constipation == 'n' ) |>
  EMP_save_var('include_df',get_result = TRUE)

MAE |>
  EMP_assay_extract('experiment') |>
  EMP_filter(primary %in% !!include_df$primary) |> 
  EMP_filter(!is.na(Metabolic_syndrome)) |>
  EMP_collapse(estimate_group='Genus',
               collapse_by = 'row') |>
  EMP_adjust_abundance(method = 'combat',
                       .factor_unwanted = 'Districts',
                       .factor_of_interest = 'Metabolic_syndrome') |>
  EMP_filter(feature_condition = !str_detect_multi(feature,'unclassified')) |>
  EMP_decostand(method = 'relative') |>
  EMP_alpha_analysis() |>
  EMP_boxplot(estimate_group = 'Metabolic_syndrome',
              select_metrics = c('shannon','simpson')) |>
  EMP_dimension_analysis(method = 'pca',scale = 'hellinger') |>
  EMP_scatterplot(show = 'p12html') |>
  EMP_identify_assay(method = 'default') |>
  EMP_diff_analysis(method = 'wilcox.test') |>
  EMP_filter(feature_condition = fdr < 0.05) |>
  EMP_save_var('microbiome_data') |>
  EMP_collapse(method = 'mean',collapse_by = 'col') |>
  EMP_heatmap_plot(scale = 'standardize',clust_row = TRUE)

(diet_data + microbiome_data + meta_data) |>
  EMP_cor_analysis() |>
  EMP_sankey_plot() |>
  EMP_network_analysis(method = 'pcor',
                       threshold = 'fdr') |>
  EMP_network_plot(show = 'node') |>
  EMP_network_plot(show = 'net',
                   shape = 'diamond', vsize = 4,
                   edge.labels = TRUE,edge.label.cex = 0.6,
                   esize = 1.5,fade = F,edge.width = 3,
                   edge.label.bg = FALSE,label.cex = 1,threshold = 0.1) 