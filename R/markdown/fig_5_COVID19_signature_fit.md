Fig 5 - COVID-19 signature fit
================

Reading data input

``` r
library(dplyr)
```

    ## Warning: replacing previous import 'vctrs::data_frame' by 'tibble::data_frame'
    ## when loading 'dplyr'

``` r
library(ggplot2)
source("../scripts/helper_functions.R")

all_contrasts <-
  readRDS(file = "../../data/mRNA_studies/all_contrasts.RDS")

all_contrasts_dict <-
  readRDS(file = "../../data/mRNA_studies/all_contrasts_dict.RDS")

filter_object <-
  readRDS("../../data/mRNA_studies/filter_object.RDS")

annotation_libraries <-
  readRDS("../../data/prior_info/gene_sets_libraries/library_list")
cell_specific_signatures <- c(annotation_libraries$iris)

COVID19_signature <- list(
  up = c("PIF1", "GUCD1", "EHD3", "TCEAL3", "BANF1"),
  down = c("ARAP2", "SLC25A46", "SLK", "ROCK2", "TVP23B", "DOCK5")
)
```

Performing signature alignment

``` r
greedy_cell_specific_alignment <-
  performa_greedy_cell_specific_alignment(
    COVID19_signature,
    cell_specific_signatures,
    all_contrasts,
    all_contrasts_dict,
    filter_object
  )
```

    ## Warning: Can't find generic `testthat_print` in package testthat to register S3 method.
    ## Can't find generic `testthat_print` in package testthat to register S3 method.
    ## Can't find generic `testthat_print` in package testthat to register S3 method.
    ## ℹ This message is only shown to developers using devtools.
    ## ℹ Do you need to update testthat to the latest version?

    ## [1] "max correlation obtained at iteration 2"

``` r
# show results
greedy_cell_specific_alignment$best_cell_signature
```

    ## [[1]]
    ## [[1]]$best_cell
    ## [1] "IRIS_PlasmaCell-FromPBMC_up"
    ## 
    ## [[1]]$best_COVID19_cell_signatures_cor2
    ## [1] 0.3411229
    ## 
    ## 
    ## [[2]]
    ## [[2]]$best_cell
    ## [1] "IRIS_PlasmaCell-FromPBMC_up+IRIS_MemoryTcell-RO-unactivated_up"
    ## 
    ## [[2]]$best_COVID19_cell_signatures_cor2
    ## [1] 0.4794752

Plotting results: scatterplot

``` r
best_alignment_AUC <-
  greedy_cell_specific_alignment$contrasts_AUC_long %>% filter(
    .id %in% c(
      "COVID19_AUC",
      "cell_signatures_AUC.IRIS_PlasmaCell-FromPBMC_up",
      "cell_signatures_AUC.IRIS_MemoryTcell-RO-unactivated_up",
      "best_cell_combo_AUC"
    )
  )

best_alignment_AUC <- best_alignment_AUC %>%
  rename(signature = .id) %>%
  dplyr::mutate(signature = gsub("COVID19_AUC", "COVID19", signature)) %>%
  dplyr::mutate(
    signature = gsub(
      "cell_signatures_AUC.IRIS_PlasmaCell-FromPBMC_up",
      "plasmablasts",
      signature
    )
  ) %>%
  dplyr::mutate(
    signature = gsub(
      "cell_signatures_AUC.IRIS_MemoryTcell-RO-unactivated_up",
      "memory_T_cells",
      signature
    )
  ) %>%
  dplyr::mutate(signature = gsub(
    "best_cell_combo_AUC",
    "plasmablasts_and_memory_T_cells",
    signature
  )) %>%
  dplyr::mutate(signature = factor(
    signature,
    levels = c(
      "COVID19",
      "plasmablasts",
      "memory_T_cells",
      "plasmablasts_and_memory_T_cells"
    )
  ))


#check statistics the AUC distributions from different signatures
lapply(split(best_alignment_AUC, 
             best_alignment_AUC$signature), 
       function(x) 
  x %>%
    group_by(class1) %>%
    summarise(perfomance = quantile(AUC, probs = c(0.25, 0.5, 0.75)), .groups = 'drop')
  )
```

    ## $COVID19
    ## # A tibble: 12 x 2
    ##    class1         perfomance
    ##    <ord>               <dbl>
    ##  1 COVID-19            0.784
    ##  2 COVID-19            0.863
    ##  3 COVID-19            0.942
    ##  4 viral               0.300
    ##  5 viral               0.393
    ##  6 viral               0.445
    ##  7 bacterial           0.173
    ##  8 bacterial           0.290
    ##  9 bacterial           0.488
    ## 10 non infectious      0.361
    ## 11 non infectious      0.488
    ## 12 non infectious      0.525
    ## 
    ## $plasmablasts
    ## # A tibble: 12 x 2
    ##    class1         perfomance
    ##    <ord>               <dbl>
    ##  1 COVID-19            0.616
    ##  2 COVID-19            0.828
    ##  3 COVID-19            0.935
    ##  4 viral               0.451
    ##  5 viral               0.527
    ##  6 viral               0.628
    ##  7 bacterial           0.317
    ##  8 bacterial           0.413
    ##  9 bacterial           0.463
    ## 10 non infectious      0.363
    ## 11 non infectious      0.489
    ## 12 non infectious      0.628
    ## 
    ## $memory_T_cells
    ## # A tibble: 12 x 2
    ##    class1         perfomance
    ##    <ord>               <dbl>
    ##  1 COVID-19           0.246 
    ##  2 COVID-19           0.369 
    ##  3 COVID-19           0.517 
    ##  4 viral              0.130 
    ##  5 viral              0.216 
    ##  6 viral              0.317 
    ##  7 bacterial          0.0355
    ##  8 bacterial          0.0986
    ##  9 bacterial          0.236 
    ## 10 non infectious     0.247 
    ## 11 non infectious     0.443 
    ## 12 non infectious     0.543 
    ## 
    ## $plasmablasts_and_memory_T_cells
    ## # A tibble: 12 x 2
    ##    class1         perfomance
    ##    <ord>               <dbl>
    ##  1 COVID-19            0.622
    ##  2 COVID-19            0.75 
    ##  3 COVID-19            0.865
    ##  4 viral               0.364
    ##  5 viral               0.427
    ##  6 viral               0.561
    ##  7 bacterial           0.138
    ##  8 bacterial           0.242
    ##  9 bacterial           0.348
    ## 10 non infectious      0.318
    ## 11 non infectious      0.477
    ## 12 non infectious      0.580

``` r
#check significance p-values for the AUC distributions from different signatures
lapply(split(best_alignment_AUC, 
             best_alignment_AUC$signature), 
       compute_AUC_distribution_p_values)
```

    ## $COVID19
    ##           class1      p_value
    ## 1       COVID-19 4.061504e-09
    ## 2          viral 7.366922e-08
    ## 3      bacterial 2.693925e-07
    ## 4 non infectious 1.840803e-02
    ## 
    ## $plasmablasts
    ##           class1      p_value
    ## 1       COVID-19 4.864470e-06
    ## 2          viral 9.735116e-01
    ## 3      bacterial 2.148082e-03
    ## 4 non infectious 3.087173e-01
    ## 
    ## $memory_T_cells
    ##           class1      p_value
    ## 1       COVID-19 9.542508e-01
    ## 2          viral 1.604151e-11
    ## 3      bacterial 5.882997e-15
    ## 4 non infectious 3.245608e-03
    ## 
    ## $plasmablasts_and_memory_T_cells
    ##           class1      p_value
    ## 1       COVID-19 2.347588e-05
    ## 2          viral 3.662863e-02
    ## 3      bacterial 3.462020e-11
    ## 4 non infectious 2.453493e-02

``` r
# scatter plot with correlation
best_alignment_AUC_scatter_df <-
  tidyr::spread(
    best_alignment_AUC %>% select(signature, study, AUC, class1, use, class3, size),
    key = signature,
    value = AUC
  )


ggplot(
  best_alignment_AUC_scatter_df,
  aes(x = plasmablasts_and_memory_T_cells, y = COVID19)
) +
  theme_Publication() +
  # geom_point(aes(color = class1), size = 4) +
  geom_point(
    size = 5,
    fill = "gray80",
    colour = "black",
    pch = 21
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "black"
  ) +
  theme(
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 12),
    legend.position = "top"
  ) +
  ylab(expression(sigma * "(COVID-19) AUCs")) +
  xlab(expression(sigma *
    "(plasmablasts + memory T) AUCs")) +
  labs(color = "contrast") +
  theme(text = element_text(size = 14))
```

![](fig_5_COVID19_signature_fit_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Plotting results: boxplots

``` r
# boxplots with performance distributions
facet_names <- as_labeller(
  c(
    COVID19 = "sigma(COVID-19)",
    plasmablasts = "sigma(plasmablasts)",
    memory_T_cells = "sigma(memory~T)",
    plasmablasts_and_memory_T_cells = "sigma(plasmablasts+memory~T)"
  ),
  default = label_parsed
)

ggplot(best_alignment_AUC, aes(x = class1, y = AUC)) +
  geom_boxplot(outlier.shape = NA, aes(color = class1), size = 1) +
  scale_color_manual(values = c("#c9daf8ff", "#fce5cdff", "#d9ead3ff", "#ead1dcff")) +
  theme_Publication() +
  #scale_fill_manual(values = c("white", paste0("gray", seq(20, 60, length.out = 3)))) +
  facet_wrap(~
  signature,
  ncol = 1,
  labeller = facet_names
  ) +
  theme(legend.position = "right") +
  xlab("contrast") +
  ylab("AUC") +
  guides(color = FALSE) +
  labs(size = "study size") +
  theme(
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  theme(strip.text.x = element_text(size = 12))
```

![](fig_5_COVID19_signature_fit_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 3.6.3 (2020-02-29)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_3.3.2 dplyr_1.0.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6      highr_0.8         pillar_1.4.4      compiler_3.6.3   
    ##  [5] plyr_1.8.6        tools_3.6.3       testthat_2.3.2    digest_0.6.25    
    ##  [9] lattice_0.20-38   nlme_3.1-144      evaluate_0.15     lifecycle_1.0.1  
    ## [13] tibble_3.0.1      gtable_0.3.0      mgcv_1.8-31       pkgconfig_2.0.3  
    ## [17] rlang_1.0.2       Matrix_1.2-18     cli_3.3.0         rstudioapi_0.11  
    ## [21] yaml_2.2.1        xfun_0.29         withr_2.2.0       stringr_1.4.0    
    ## [25] knitr_1.39        generics_0.0.2    vctrs_0.4.1       grid_3.6.3       
    ## [29] tidyselect_1.1.0  glue_1.6.2        R6_2.4.1          fansi_0.4.1      
    ## [33] rmarkdown_2.3     farver_2.0.3      purrr_0.3.4       tidyr_1.1.0      
    ## [37] magrittr_2.0.1    splines_3.6.3     scales_1.1.1      ellipsis_0.3.2   
    ## [41] htmltools_0.5.1.1 ggthemes_4.2.0    colorspace_1.4-1  labeling_0.3     
    ## [45] utf8_1.1.4        stringi_1.4.6     munsell_0.5.0     crayon_1.3.4
