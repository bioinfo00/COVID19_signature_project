Fig 3 and Fig S6 - COVID-19 signature validation
================

This Mardown file generates the figure panels related to the validation
of the COVID-19 specific signature.

Panel A: generating ROC curves when applying the signature to GSE163151.

``` r
library(dplyr)
library(ggplot2)
source("../scripts/helper_functions.R")

COVID19_signature_up <- c("PIF1", "GUCD1", "EHD3", "TCEAL3", "BANF1")
COVID19_signature_down <- c("ARAP2", "SLC25A46", "SLK", "ROCK2", "TVP23B", "DOCK5")

# analysis of GSE163151 which has COVID-19, viral, and bacterial samples
dataset <- readRDS(file = "../../data/mRNA_studies/GSE163151_preprocessed.RDS")
COVID19_signature_up_expression <- dataset$X[COVID19_signature_up, ]
COVID19_signature_down_expression <- dataset$X[COVID19_signature_down, ]

sample_score <- scale(
  apply(COVID19_signature_up_expression, 2, function(x) {
    geom_mean(x)
  }) -
    apply(COVID19_signature_down_expression, 2, function(x) {
      geom_mean(x)
    })
)
```

    ## Error in get(genname, envir = envir) : object 'testthat_print' not found

``` r
dataset$y[, "disease state:ch1"] <- factor(
  dataset$y[, "disease state:ch1"],
  levels = c(
    "Donor control",
    "COVID-19",
    "Viral acute respiratory illness",
    "Bacterial sepsis"
  )
)


score_df <- data.frame(
  score = sample_score,
  group = dataset$y[, "disease state:ch1"]
) %>%
  dplyr::mutate(group = relevel(group, ref = "Donor control"))


# ROC curve related to COVID-19 vs control
covid_vs_control_roc_curve <- get_roc_curve(score_df,
  control_group = "Donor control",
  positive_group = "COVID-19"
)

# ROC curve related to viral vs control
viral_vs_control_roc_curve <- get_roc_curve(score_df,
  control_group = "Donor control",
  positive_group = "Viral acute respiratory illness"
)

# ROC curve related to bacterial vs control
bact_vs_control_roc_curve <- get_roc_curve(score_df,
  control_group = "Donor control",
  positive_group = "Bacterial sepsis"
)

p_roc <- pROC::ggroc(
  list(
    COVID = covid_vs_control_roc_curve,
    viral = viral_vs_control_roc_curve,
    bacterial = bact_vs_control_roc_curve
  ),
  size = 2
) + theme_minimal() +
  scale_colour_manual(values = c("#c9daf8ff", "#fce5cdff", "#d9ead3ff")) +
  annotate(
    "text",
    x = 0.72,
    y = 0.9,
    label = paste(
      "COVID-19 contrast \n AUC=",
      round(pROC::auc(covid_vs_control_roc_curve), 2)
    ),
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.52,
    y = 0.45,
    label = paste(
      "viral respiratory contrast \n AUC=",
      round(pROC::auc(viral_vs_control_roc_curve), 2)
    ),
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.4,
    y = 0.1,
    label = paste(
      "bact
erial sepsis contrast \n AUC=",
      round(pROC::auc(bact_vs_control_roc_curve), 3)
    ),
    hjust = 0
  ) +
  theme_Publication() + theme(legend.position = "none") +
  ggtitle("GSE163151 validation") + theme(text = element_text(size = rel(4.5))) +
  theme(plot.title = element_text(size = rel(3.5)))

p_roc
```

![](fig_3_and_fig_S6_signature_validation_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Panel B: generating ROC curves when applying the signature to GSE149689

``` r
dataset <- readRDS(file = "../../data/mRNA_studies/COVID19_validation/GSE149689_pseudobulk_MI.RDS")
COVID19_signature_up_expression <- dataset$expr[COVID19_signature_up, ]
COVID19_signature_down_expression <- dataset$expr[COVID19_signature_down, ]

sample_score <- scale(
  apply(COVID19_signature_up_expression, 2, function(x) {
    geom_mean(x)
  }) -
    apply(COVID19_signature_down_expression, 2, function(x) {
      geom_mean(x)
    })
)


score_df <- data.frame(
  score = sample_score,
  group = dataset$pheno$covid_status
) %>%
  dplyr::mutate(group = relevel(group, ref = "Healthy"))


# ROC curve related to COVID-19 vs control
covid_vs_control_roc_curve <- get_roc_curve(score_df,
                                            control_group = "Healthy",
                                            positive_group = "COVID-19"
)

# ROC curve related to viral vs control
viral_vs_control_roc_curve <- get_roc_curve(score_df,
                                            control_group = "Healthy",
                                            positive_group = "Flu"
)


p_roc_GSE149689 <- pROC::ggroc(
  list(
    COVID = covid_vs_control_roc_curve,
    viral = viral_vs_control_roc_curve
  ),
  size = 2
) + theme_minimal() +
  scale_colour_manual(values = c("#c9daf8ff", "#fce5cdff")) +
  annotate(
    "text",
    x = 0.72,
    y = 0.9,
    label = paste(
      "COVID-19 contrast \n AUC=",
      round(pROC::auc(covid_vs_control_roc_curve), 2)
    ),
    hjust = 0
  ) +
  annotate(
    "text",
    x = 0.46,
    y = 0.45,
    label = paste(
      "viral respiratory contrast \n AUC=",
      round(pROC::auc(viral_vs_control_roc_curve), 2)
    ),
    hjust = 0
  ) +
  theme_Publication() + theme(legend.position = "none") +
  ggtitle("GSE149689 validation") + theme(text = element_text(size = rel(4.5))) +
  theme(plot.title = element_text(size = rel(3.5)))
p_roc_GSE149689
```

![](fig_3_and_fig_S6_signature_validation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Panel C: AUC distribution including all validation studies

``` r
filter_object <-
  readRDS("../../data/mRNA_studies/filter_object.RDS")

all_contrasts <-
  readRDS(file = "../../data/mRNA_studies/all_contrasts.RDS")
all_contrasts_dict <-
  readRDS(file = "../../data/mRNA_studies/all_contrasts_dict.RDS")

COVID19_signature = list(up = COVID19_signature_up, 
                         down = COVID19_signature_down)

COVID19_validation_contrasts_AUC <- compute_contrasts_AUCs(
  COVID19_signature,
  filter_object,
  all_contrasts[which(all_contrasts_dict$use == 'validation')],
  all_contrasts_dict %>% filter(use == 'validation')
)


#generating summary statistics (first quartile, median, third quartile)
COVID19_validation_contrasts_AUC %>%
  group_by(class1, use) %>%
  summarise(perfomance = quantile(AUC, probs = c(0.25, 0.5, 0.75)))
```

    ## # A tibble: 12 x 3
    ## # Groups:   class1, use [4]
    ##    class1         use        perfomance
    ##    <ord>          <fct>           <dbl>
    ##  1 COVID-19       validation      0.757
    ##  2 COVID-19       validation      0.805
    ##  3 COVID-19       validation      0.919
    ##  4 viral          validation      0.258
    ##  5 viral          validation      0.356
    ##  6 viral          validation      0.485
    ##  7 bacterial      validation      0.182
    ##  8 bacterial      validation      0.273
    ##  9 bacterial      validation      0.519
    ## 10 non infectious validation      0.340
    ## 11 non infectious validation      0.478
    ## 12 non infectious validation      0.509

``` r
p_boxplot = generate_AUC_boxplot(COVID19_validation_contrasts_AUC)

p_boxplot
```

![](fig_3_and_fig_S6_signature_validation_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Fig S6: AUC distributions produced by previously published COVID-19
signatures

``` r
published_signatures = readRDS('../../data/published_signatures/published_signatures.RDS')


published_signatures_validation_AUC = lapply(published_signatures, 
                                             function(x) compute_contrasts_AUCs(x, filter_object,
                                          all_contrasts[which(all_contrasts_dict$use == 'validation')],
                                               all_contrasts_dict %>% filter(use == 'validation')
                                             ))

#from Thair, exclude the discovery study
published_signatures_validation_AUC$Thair = published_signatures_validation_AUC$Thair %>% 
  filter(!grepl('GSE152641', study))


#generating summary statistics (first quartile, median, third quartile)
lapply(published_signatures_validation_AUC, function(x) 
  x %>%
    group_by(class1, use) %>%
    summarise(perfomance = quantile(AUC, probs = c(0.25, 0.5, 0.75)), .groups = 'drop')
  %>% ungroup()
  )
```

    ## $Thair
    ## # A tibble: 12 x 3
    ##    class1         use        perfomance
    ##    <ord>          <fct>           <dbl>
    ##  1 COVID-19       validation      0.592
    ##  2 COVID-19       validation      0.833
    ##  3 COVID-19       validation      0.979
    ##  4 viral          validation      0.480
    ##  5 viral          validation      0.580
    ##  6 viral          validation      0.784
    ##  7 bacterial      validation      0.647
    ##  8 bacterial      validation      0.699
    ##  9 bacterial      validation      0.75 
    ## 10 non infectious validation      0.452
    ## 11 non infectious validation      0.528
    ## 12 non infectious validation      0.596
    ## 
    ## $Lee
    ## # A tibble: 12 x 3
    ##    class1         use        perfomance
    ##    <ord>          <fct>           <dbl>
    ##  1 COVID-19       validation      0.694
    ##  2 COVID-19       validation      0.861
    ##  3 COVID-19       validation      0.929
    ##  4 viral          validation      0.443
    ##  5 viral          validation      0.596
    ##  6 viral          validation      0.741
    ##  7 bacterial      validation      0.694
    ##  8 bacterial      validation      0.843
    ##  9 bacterial      validation      0.888
    ## 10 non infectious validation      0.553
    ## 11 non infectious validation      0.598
    ## 12 non infectious validation      0.669
    ## 
    ## $McClain
    ## # A tibble: 12 x 3
    ##    class1         use        perfomance
    ##    <ord>          <fct>           <dbl>
    ##  1 COVID-19       validation      0.571
    ##  2 COVID-19       validation      0.803
    ##  3 COVID-19       validation      0.838
    ##  4 viral          validation      0.558
    ##  5 viral          validation      0.827
    ##  6 viral          validation      0.894
    ##  7 bacterial      validation      0.271
    ##  8 bacterial      validation      0.453
    ##  9 bacterial      validation      0.541
    ## 10 non infectious validation      0.388
    ## 11 non infectious validation      0.457
    ## 12 non infectious validation      0.515
    ## 
    ## $Aschenbrenner
    ## # A tibble: 12 x 3
    ##    class1         use        perfomance
    ##    <ord>          <fct>           <dbl>
    ##  1 COVID-19       validation      0.777
    ##  2 COVID-19       validation      0.922
    ##  3 COVID-19       validation      0.986
    ##  4 viral          validation      0.662
    ##  5 viral          validation      0.818
    ##  6 viral          validation      0.875
    ##  7 bacterial      validation      0.734
    ##  8 bacterial      validation      0.960
    ##  9 bacterial      validation      0.983
    ## 10 non infectious validation      0.491
    ## 11 non infectious validation      0.631
    ## 12 non infectious validation      0.796

``` r
published_signatures_validation_AUC_df = plyr::ldply(published_signatures_validation_AUC, 
                                                     rbind)
published_signatures_validation_AUC_boxplot = generate_AUC_boxplot(published_signatures_validation_AUC_df) +
  facet_wrap(~.id, ncol = 1)

published_signatures_validation_AUC_boxplot
```

![](fig_3_and_fig_S6_signature_validation_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

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
    ##  [1] Rcpp_1.0.4.6      pillar_1.4.4      compiler_3.6.3    plyr_1.8.6       
    ##  [5] tools_3.6.3       testthat_2.3.2    digest_0.6.25     evaluate_0.14    
    ##  [9] lifecycle_0.2.0   tibble_3.0.1      gtable_0.3.0      pkgconfig_2.0.3  
    ## [13] rlang_0.4.11      cli_2.0.2         yaml_2.2.1        xfun_0.15        
    ## [17] withr_2.2.0       stringr_1.4.0     knitr_1.29        generics_0.0.2   
    ## [21] pROC_1.16.2       vctrs_0.3.1       grid_3.6.3        tidyselect_1.1.0 
    ## [25] glue_1.4.1        R6_2.4.1          fansi_0.4.1       rmarkdown_2.3    
    ## [29] purrr_0.3.4       farver_2.0.3      magrittr_2.0.1    scales_1.1.1     
    ## [33] ellipsis_0.3.1    htmltools_0.5.1.1 ggthemes_4.2.0    assertthat_0.2.1 
    ## [37] colorspace_1.4-1  labeling_0.3      utf8_1.1.4        stringi_1.4.6    
    ## [41] munsell_0.5.0     crayon_1.3.4
