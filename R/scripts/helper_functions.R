theme_Publication <-
  function(base_size = 14,
           base_family = "Helvetica") {
    # ggplot graphic options for pretty graphics
    (
      ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
      + theme(
          plot.title = element_text(
            face = "bold",
            size = rel(1.2),
            hjust = 0.5
          ),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(size = rel(0.9)),
          axis.title.y = element_text(angle = 90, vjust = 2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour = "#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "left",
          legend.key.size = unit(0.2, "cm"),
          legend.spacing = unit(0, "cm"),
          plot.margin = unit(c(10, 5, 5, 5), "mm"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
          strip.text = element_text(face = "bold")
        )
    )
  }

create_signature_list <- function(signature_gene_char_vector) {
  # turns a signature represented as a character vector into a list
  testthat::expect_true(
    class(signature_gene_char_vector) == "character",
    label = "signature vector is a character vector",
    info = "input should be a character vector"
  )
  testthat::expect_false(any(is.na(signature_gene_char_vector)),
    label = "signature vector has NAs",
    info = "input should not contain NAs"
  )
  signature_list <- list(
    up = signature_gene_char_vector,
    down = ""
  )
  return(signature_list)
}

reverse_signature_list <- function(signature_list) {
  # reverses a signature list by swapping its up- and down-regulated genes
  reversed_signature_list <- list(
    up = signature_list$down,
    down = signature_list$up
  )
  return(reversed_signature_list)
}

create_cell_type_signature_list <- function(cell_type_annotation) {
  # creates a list of cell type signature lists from an cell type annotation library
  # it includes signature lists corresponding to up- and down-regulation of each cell type
  # in the annotation library
  cell_type_signature_lists <-
    lapply(cell_type_annotation, create_signature_list)
  names(cell_type_signature_lists) <-
    paste(names(cell_type_signature_lists),
      "up",
      sep = "_"
    )

  reversed_cell_type_signature_lists <-
    lapply(cell_type_signature_lists, reverse_signature_list)
  names(reversed_cell_type_signature_lists) <-
    paste(names(reversed_cell_type_signature_lists),
      "down",
      sep = "_"
    )

  cell_type_signature_lists <- c(
    cell_type_signature_lists,
    reversed_cell_type_signature_lists
  )

  return(cell_type_signature_lists)
}

combine_signature_lists <-
  function(signature_list_1, signature_list_2) {
    # creates a signature list for the combination of two signature lists
    # from two cell types
    up_union <- union(signature_list_1$up, signature_list_2$up)
    down_union <-
      union(signature_list_1$down, signature_list_2$down)
    up_down_intersection <- intersect(up_union, down_union)

    # remove from the up- and down-regulated components the discordant genes
    up_union <- setdiff(up_union, up_down_intersection)
    down_union <- setdiff(down_union, up_down_intersection)

    # expect_false(length(up_union) == 0 && length(down_union) == 0,
    #   label = "distinct up- and down-regulated genes"
    # )

    # check if the same genes are both up- and down-regulated
    combined_signature <- list(
      up = up_union,
      down = down_union
    )

    return(combined_signature)
  }

geom_mean <- function(x) {
  # computes the geometric mean of a numeric vector x
  testthat::expect_true(class(x) == "numeric", info = "input should be numeric")

  NA_values <- sum(is.na(x))
  if (sum(NA_values) > 0) {
    warning("one or more input values are missing")
  }

  return(exp(mean(log(x), na.rm = T)))
}

compute_sample_scores <- function(datasetObject, genes) {
  # computes the geometric mean of a set of genes in a
  # datasetObject consistent with the meta-integrator package
  genes_idx <- which(row.names(datasetObject$expr) %in% genes)
  if (length(genes_idx) == 0) {
    return(rep(0, ncol(datasetObject$expr)))
  }

  gene_expr <- datasetObject$expr[genes_idx, ]
  # in case only one gene is present, take the vector itself
  if (is.null(nrow(gene_expr))) {
    sample_scores <- gene_expr
  }
  if (!is.null(nrow(gene_expr))) {
    sample_scores <- apply(gene_expr, 2, geom_mean)
  }

  return(sample_scores)
}

calculate_score_fast <- function(filter_object,
                                 datasetObject,
                                 suppressMessages = T) {
  # this function compute the score of each sample in a study with respect
  # to a signature. The score is defined as in the meta-integrator package:
  # geometric mean of the up-regulated genes - geometric mean of the
  # down-regulated genes in the signature

  # pull out positive (up-regulated) and negative (down_regulated)
  # genes from signature
  pos.genes <- filter_object$posGeneNames
  neg.genes <- filter_object$negGeneNames

  # intersect these with the genes present in the expression matrix
  genes_present <- row.names(datasetObject$expr)
  pos.genes_present <- intersect(pos.genes, genes_present)
  neg.genes_present <- intersect(neg.genes, genes_present)

  posScore <-
    compute_sample_scores(datasetObject, pos.genes_present)
  negScore <-
    compute_sample_scores(datasetObject, neg.genes_present)

  totalScore <- posScore - negScore

  if (sum(abs(totalScore)) != 0) {
    totalScore <- as.numeric(scale(totalScore))
  }

  if (!suppressMessages) {
    cat(
      "Used ",
      length(pos.genes_present),
      "of ",
      length(pos.genes),
      " pos genes, and ",
      length(neg.genes_present),
      " of ",
      length(neg.genes),
      " neg genes \n"
    )
  }
  return(totalScore)
}

auroc <- function(bool, score) {
  # fast implementation of ROC AUC calculation
  # bool is a vector of binay labels
  # score is a vector of scores

  n_classes <- length(unique(bool))
  # response vector must have 2 levels or setting AUC to NA
  if (n_classes != 2) {
    print("response must have 2 levels: setting AUC to NA")
    return(NA)
  }
  # score vecor must have >=2 levels or setting AUC to NA
  n_score_values <- length(unique(score))
  if (n_score_values < 2) {
    return(NA)
  }

  n1 <- sum(!bool)
  n2 <- sum(bool)
  U <- sum(rank(score)[!bool], na.rm = T) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

compute_contrasts_AUCs <- function(signature,
                                   filter_object,
                                   contrasts,
                                   contrasts_dict) {
  filter_object_temp <- filter_object
  filter_object_temp$posGeneNames <- signature$up
  filter_object_temp$negGeneNames <- signature$down

  scores <- lapply(contrasts, function(x) {
    data.frame(
      score = calculate_score_fast(filter_object_temp, x),
      group = x$class
    )
  })

  aucs <- sapply(scores, function(x) {
    auroc(x$group, x$score)
  })

  contrasts_dict_with_AUCs <- contrasts_dict

  contrasts_dict_with_AUCs$AUC <- aucs

  return(contrasts_dict_with_AUCs)
}


performa_greedy_cell_specific_alignment <-
  function(COVID19_signature,
           cell_specific_signatures,
           all_contrasts,
           all_contrasts_dict,
           filter_object,
           max_iter = 5) {
    COVID19_AUCs <- compute_contrasts_AUCs(
      COVID19_signature,
      filter_object,
      contrasts = all_contrasts,
      contrasts_dict = all_contrasts_dict
    )

    cell_specific_signatures_lists <-
      create_cell_type_signature_list(cell_specific_signatures)

    # initialization
    best_cell_signature <- list()
    best_COVID19_cell_signatures_cor2 <- 0
    temp_cell_signatures <- cell_specific_signatures_lists

    for (k in 1:max_iter) {
      temp_cell_signatures_AUC <- lapply(
        temp_cell_signatures,
        function(x) {
          compute_contrasts_AUCs(x,
            filter_object,
            contrasts = all_contrasts,
            contrasts_dict = all_contrasts_dict
          )
        }
      )

      # get AUCs for individual cell types
      if (k == 1) {
        cell_signatures_AUC <- temp_cell_signatures_AUC
      }

      COVID19_cell_signatures_cor2 <-
        sapply(
          temp_cell_signatures_AUC,
          function(x) {
            cor(x$AUC, COVID19_AUCs$AUC)^2
          }
        )
      current_COVID19_cell_signatures_cor2 <-
        max(COVID19_cell_signatures_cor2, na.rm = T)

      if (current_COVID19_cell_signatures_cor2 > best_COVID19_cell_signatures_cor2) {
        best_COVID19_cell_signatures_cor2 <-
          current_COVID19_cell_signatures_cor2
        best_cell <- names(which.max(COVID19_cell_signatures_cor2))
        best_cell_signature[[k]] <- list(
          best_cell = best_cell,
          best_COVID19_cell_signatures_cor2 = best_COVID19_cell_signatures_cor2
        )

        temp_cell_signatures <- lapply(
          temp_cell_signatures,
          function(x) {
            combine_signature_lists(temp_cell_signatures[[best_cell]], x)
          }
        )
        names(temp_cell_signatures) <-
          paste(best_cell, names(temp_cell_signatures), sep = "+")
        best_cell_combo_AUC <- temp_cell_signatures_AUC[[best_cell]]
      } else {
        print(paste("max correlation obtained at iteration", k - 1))

        break
      }
    }

    contrasts_AUC_long <- c(
      list(COVID19_AUC = COVID19_AUCs),
      cell_signatures_AUC = cell_signatures_AUC,
      list(best_cell_combo_AUC = best_cell_combo_AUC)
    ) %>% plyr::ldply(rbind)

    return(
      list(
        best_cell_signature = best_cell_signature,
        best_COVID19_cell_signatures_cor2 = best_COVID19_cell_signatures_cor2,
        contrasts_AUC_long = contrasts_AUC_long
      )
    )
  }


get_roc_curve <- function(score_df, control_group, positive_group) {
  positive_vs_control <-
    score_df %>% filter(group %in% c(control_group, positive_group))
  covid_control_roc <-
    pROC::roc(droplevels(positive_vs_control$group),
      positive_vs_control$score,
      direction = "<"
    )
  return(covid_control_roc)
}



plot_COVID19_signature_scores <- function(study_accession,
                                          validation_study_files,
                                          sex_field,
                                          age_field,
                                          severity_field) {
  validation_study_files <- list.files(path = validation_study_path)

  # get expression matrix
  validation_study_file_idx <-
    grep(study_accession, validation_study_files)
  study_file <- validation_study_files[validation_study_file_idx]
  study <-
    readRDS(file = paste0(validation_study_path, "/", study_file))
  expression_matrix <- study[[1]]

  COVID19_signature_up_expression <-
    expression_matrix[COVID19_signature_up, ]
  COVID19_signature_down_expression <-
    expression_matrix[COVID19_signature_down, ]

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
    group = study[[2]]
  )

  boxplots <- list(
    sex = NULL,
    age = NULL,
    severity = NULL
  )

  if (!is.na(sex_field)) {
    boxplots$sex <- ggplot(
      score_df %>% tidyr::drop_na(sex_field),
      aes_string(x = sex_field, y = "score")
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter() +
      ggtitle(study_accession) +
      theme_Publication() +
      xlab("sex") +
      ylab("") +
      theme(text = element_text(size = 16))
  }

  if (!is.na(age_field)) {
    score_df[, age_field] <- as.numeric(score_df[, age_field])
    score_df[, age_field] <-
      cut(score_df[, age_field], breaks = seq(0, 100, by = 20))
    boxplots$age <- ggplot(
      score_df %>% tidyr::drop_na(age_field),
      aes_string(x = age_field, y = "score")
    ) +
      geom_boxplot() +
      geom_jitter() +
      ggtitle(study_accession) +
      theme_Publication() +
      xlab("age") +
      ylab("") +
      theme(text = element_text(size = 16))
  }

  if (!is.na(severity_field)) {
    boxplots$severity <- ggplot(
      score_df %>% tidyr::drop_na(severity_field),
      aes_string(x = severity_field, y = "score")
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter() +
      ggtitle(study_accession) +
      theme_Publication() +
      xlab("severity") +
      ylab("") +
      theme(text = element_text(size = 20))
  }

  return(boxplots)
}


# not used in markdown files
# extract_population_from_GA_object <- function(GA_object,
#                                               min_fitness = 0.0) {
#   GA_population <- data.frame(GA_object@population,
#                               fitness = as.numeric(GA_object@fitness))
#
#   GA_population <- GA_population %>%
#     dplyr::filter(fitness > min_fitness) %>%
#     dplyr::arrange(desc(fitness))
#
#   GA_population <- GA_population %>%
#     dplyr::select(-fitness) %>%
#     as.matrix()
#
#   return(GA_population)
# }


get_signature_from_binary_solution <- function(binary_solution, filter_object) {
  pool_of_genes <- names(binary_solution)
  signature_genes <- pool_of_genes[which(binary_solution == 1)]
  up <- intersect(signature_genes, filter_object$posGeneNames)
  down <- intersect(signature_genes, filter_object$negGeneNames)
  signature <- list(
    up = up,
    down = down
  )
  return(signature)
}




compute_multi_objective_fitness <- function(contrasts_AUC_df,
                                            AUC_field,
                                            class_field,
                                            detection_pattern = "COVID") {
  # takes a dataframe contrasts_AUC_df that has a set of contrasts for detection
  # and cross-reactivity
  # the dataframe should contain two fields: 1) an AUC produced by a given signature;
  # 2) a factor classifying the contrasts
  # AUC = sym(AUC_field)
  # class = sym(class_field)


  detection_fitness_components <- compute_fitness_components(
    contrasts_AUC_df,
    AUC_field = AUC_field,
    class_field = class_field,
    detection_pattern = detection_pattern,
    type = "detection"
  )

  cross_reactivity_fitness_components <- compute_fitness_components(
    contrasts_AUC_df,
    AUC_field = AUC_field,
    class_field = class_field,
    detection_pattern = detection_pattern,
    type = "cross_reactivity"
  )

  multi_objective_fitness <- rbind(
    detection_fitness_components,
    cross_reactivity_fitness_components
  ) %>% data.frame()


  return(multi_objective_fitness)
}


compute_fitness_components <- function(contrasts_AUC_df,
                                       AUC_field,
                                       class_field,
                                       detection_pattern = "COVID",
                                       type) {
  
  if (type == "detection") {
    fitness_components <- contrasts_AUC_df %>%
      dplyr::filter(grepl(detection_pattern, !!sym(class_field))) %>%
      dplyr::group_by(across(class_field)) %>%
      dplyr::summarise(
        fitness = compute_fitness(AUC, type = type),
        .groups = "drop"
      )
  } else {
    fitness_components <- contrasts_AUC_df %>%
      dplyr::filter(!grepl(detection_pattern, !!sym(class_field))) %>%
      dplyr::group_by(across(class_field)) %>%
      dplyr::summarise(
        fitness = compute_fitness(AUC, type = type),
        .groups = "drop"
      )
  }

  return(fitness_components)
}

compute_fitness <- function(AUC_vector, type) {
  # computes fitness by aggregating the AUCs of contrasts within the same class
  # the aggregation works differently for detection and cross-reactivity contrasts
  fitness <- ifelse(type == "detection",
    min(AUC_vector, na.rm = T),
    1 - 2 * max(ifelse(AUC_vector > 0.5, AUC_vector - 0.5, 0), na.rm = T)
  )

  return(fitness)
}


compute_mean_distance_from_utopia <- function(fitness_vector) {
  # given a fitness_vector, it computes the distance from
  # the utopia point (all fitness components = 1)

  # removing NAs from vector, if any
  fitness_vector <- na.omit(fitness_vector)
  n_of_fitness_components <- length(fitness_vector)

  utopia <- rep(1, n_of_fitness_components)
  distance_from_utopia <- sqrt(mean(fitness_vector - utopia)^2)

  return(distance_from_utopia)
}


get_solution_stability <- function(binary_solution_matrix) {
  gene_stability <- apply(binary_solution_matrix, 2, function(x) {
    sum(x) / length(x)
  }) %>%
    reshape2::melt(gene_stability) %>%
    tibble::rownames_to_column("gene") %>%
    rename(stability = value)

  solution_stability <- apply(
    binary_solution_matrix, 1,
    function(x) {
      mean(gene_stability[
        which(x == 1),
        "stability"
      ])
    }
  )

  return(solution_stability)
}

get_signature_mean_gene_overlap <- function(candidate_signatures,
                                            contrasts) {
  mean_gene_overlap <- vector(length = length(candidate_signatures))
  for (k in 1:length(candidate_signatures)) {
    signature_genes <- unlist(candidate_signatures[[k]])
    mean_gene_overlap[k] <- mean(sapply(
      contrasts,
      function(x) {
        length(intersect(
          x$keys, signature_genes
        ))
      }
    ))
  }
  return(mean_gene_overlap)
}

get_distances_from_utopia <- function(signatures,
                                      filter_object,
                                      contrasts,
                                      contrasts_dict) {
  contrasts_AUC <- lapply(signatures, function(x) {
    compute_contrasts_AUCs(
      x,
      filter_object,
      contrasts,
      contrasts_dict
    )
  })

  signature_fitness <- lapply(contrasts_AUC, function(x) {
    compute_multi_objective_fitness(x,
      class_field = "class2"
    )
  })


  dist_from_utopia <- sapply(
    signature_fitness,
    function(x) {
      compute_mean_distance_from_utopia(x$fitness)
    }
  )

  return(dist_from_utopia)
}


get_AUC_distribution_boxplot <- function(contrasts_AUC_df) {
  # generate a boxplot with distribution of AUC corresponding to the
  # different study classes (COVID-19, viral, bacterial, non-infectious)
  contrasts_AUC_df$respiratory <- "not"
  contrasts_AUC_df$respiratory[grep("Resp", contrasts_AUC_df$class)] <-
    "yes"

  contrasts_AUC_df$respiratory[grep(
    "COVID",
    contrasts_AUC_df$class
  )] <- "yes"

  p_boxplot <- ggplot(
    contrasts_AUC_df,
    aes(x = class1, y = AUC)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(
      size  = size,
      color = class1,
      shape = respiratory
    )) +
    scale_color_manual(values = c("#c9daf8ff", "#fce5cdff", "#d9ead3ff", "#ead1dcff")) +
    theme_Publication() +
    theme(legend.position = "right") +
    xlab("contrast") +
    ylab("AUC") +
    guides(color = FALSE) +
    labs(size = "study size") +
    theme(
      legend.key.size = unit(2, "line"),
      legend.text = element_text(size = 16),
      legend.position = "top"
    ) +
    theme(strip.text.x = element_text(size = 12)) +
    guides(
      size = guide_legend(order = 2),
      shape = guide_legend(order = 1, override.aes = list(size = 5))
    ) +
    guides(guide_legend(nrow = 2, byrow = TRUE)) +
    theme(text = element_text(size = 20)) +
    theme(legend.direction = "horizontal", legend.box = "vertical") +
    annotate(
      "text",
      x = 1,
      y = 1.1,
      label = "detection",
      size = 6
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    annotate(
      "text",
      x = 3,
      y = 1.1,
      label = "cross-reactivity",
      size = 6
    ) +
    geom_segment(aes(
      x = 0.6,
      y = 1.05,
      xend = 1.4,
      yend = 1.05
    )) +
    geom_segment(aes(
      x = 1.6,
      y = 1.05,
      xend = 4.4,
      yend = 1.05
    ))

  return(p_boxplot)
}


compute_AUC_distribution_p_values = function(contrasts_AUC){
  # compute significance p-values based on hypothesis testing
  
  # get p-values for COVID-19 detection class
  # in this case, the null hypothesis is mean AUC <= 0.5
  detection_pvalue = contrasts_AUC %>% 
    dplyr::filter(class1 == 'COVID-19') %>% dplyr::group_by(class1) %>% 
    summarise(p_value = t.test(AUC, mu = 0.5, alternative = 'greater')$p.value, 
              .groups = 'drop')
  
  # get p-values for cross-reactivity in the non-COVID-19 classes
  # in this case, the null hypothesis is mean AUC >= 0.5
  cross_reactivity_pvalues = contrasts_AUC %>% 
    dplyr::filter(!class1 == 'COVID-19') %>% 
    dplyr::group_by(class1) %>% 
    summarise(p_value = t.test(AUC, mu = 0.5, alternative = 'less')$p.value, 
              .groups = 'drop')
  
  combined_p_values = rbind(detection_pvalue, cross_reactivity_pvalues) %>% 
    data.frame()
  
  return(combined_p_values)
  
}
