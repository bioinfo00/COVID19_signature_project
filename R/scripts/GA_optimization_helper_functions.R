initialize_population = function(population_size, solution_size, 
                                 expected_signature_size, filter_object){
  #generate an initial population of solutions for GA
  #make sure that the initial population consists of feasible solutions:
  #1) right size (number of 1's) and 2) balance between up- and down-regulated genes
  n_up = length(filter_object$posGeneNames)
  n_down = length(filter_object$negGeneNames)
  initial_population = matrix(0, nrow = population_size, ncol = solution_size)
  
  for (k in 1:population_size){
    
    sample_up_idx = sample(1:n_up, expected_signature_size/2)
    sample_down_idx = sample((n_up + 1):solution_size, expected_signature_size/2)
    initial_population[k, c(sample_up_idx, sample_down_idx)] = 1
    
  }
  
  return(initial_population)
}

select_initial_population_with_prior_info = function(initial_population, 
                                                     prior_info_matrix, 
                                                     top_n = 100){
  #get a randomly generated population of solutions and select the top_n
  #solutions with the largest mean projection along the prior info components
  initial_population_projection = prior_info_matrix %*% t(initial_population)
  mean_projection = apply(initial_population_projection, 2, mean)
  
  population_with_best_projection_idx = order(mean_projection, decreasing = T)[1:top_n]
  initial_population_with_prior_info = initial_population[population_with_best_projection_idx, ]
  
  return(initial_population_with_prior_info)
  
}


create_weight_grid = function(from = 0, to = 1, by = 0.5){
  #generate convex combinations of the 6 weights
  
  weight_grid = data.frame(expand.grid(seq(from = from, to = to, by = 1),
                                       seq(from = from, to = to, by = 1),
                                       seq(from = from, to = to, by = by),
                                       seq(from = from, to = to, by = by),
                                       seq(from = from, to = to, by = by),
                                       seq(from = from, to = to, by = 1)))
  
  #remove the point with all zero's
  weight_grid = weight_grid[-which(apply(weight_grid, 1, sum) == 0), ]
  #remove the point with all zero's but prior info
  weight_grid = weight_grid[-which(apply(weight_grid[, 1:5], 1, sum) == 0), ]
  
  weight_grid = t(apply(weight_grid, 1, function(x) x/sum(x)))
  weight_grid = weight_grid[!duplicated(weight_grid), ]
  
  row.names(weight_grid) = apply(weight_grid, 1, function(x) paste0('x_', paste(round(x, 2), collapse = '_')))
  colnames(weight_grid) = c('weight_COVID_vs_healthy', 
                            'weight_COVID_vs_infection', 
                            'weight_Non.Microbial',
                            'weight_Other', 
                            'weight_Resp',
                            'weight_prior_info')
  return(weight_grid)
}


compute_fitness_from_binary = function(binary_vector,
                                       pool_of_genes,
                                       filter_object,
                                       weight_COVID_vs_healthy, 
                                       weight_COVID_vs_infection, 
                                       weight_Non.Microbial,
                                       weight_Other, 
                                       weight_Resp,
                                       prior_info_matrix,
                                       weight_prior_info,
                                       max_signature_size){
  #given a signature encoded as a binary vector, this function
  #evaluates its fitness based on COVID and non COVID
  #contrasts and on the amount of prior info. 
  #Different classes of contrasts can be weighted differently by modifying the weights
  
  binary_vector_genes = pool_of_genes[which(binary_vector == 1)]
  
  temp_filter_object = filter_object
  temp_filter_object$posGeneNames = intersect(filter_object$posGeneNames, binary_vector_genes)
  temp_filter_object$negGeneNames = intersect(filter_object$negGeneNames, binary_vector_genes)
  
  n_positive = length(temp_filter_object$posGeneNames)
  n_negative = length(temp_filter_object$negGeneNames)
  
  #make sure that the signature satisfy constraints on imbalance and size
  is_signature_ok = check_signature(n_positive, n_negative, max_signature_size)
  if (!is_signature_ok) return(-2)
  
  multi_objective_vector = compute_multi_objective_vector(filter_object = temp_filter_object)
  
  total_fitness = scalarize_multi_objective_vector(multi_objective_vector, 
                                                   weight_COVID_vs_healthy, 
                                                   weight_COVID_vs_infection, 
                                                   weight_Non.Microbial, 
                                                   weight_Other, 
                                                   weight_Resp)
  
  #add term measureing the amount of prior info contained in the signature  
  if (!is.na(prior_info_matrix)[1]){
    fitness_prior_info = mean(prior_info_matrix %*% binary_vector/(sqrt(sum(binary_vector^2))))
    total_fitness = total_fitness + weight_prior_info*fitness_prior_info
    
  }
  
  return(total_fitness)
}


calculate_score_fast = function (filter_object, datasetObject, suppressMessages = T){
  #this function compute the score of each sample in a study with respect
  #to a signature. The score is defined as in the meta-integrator package:
  #geometric mean of the up-regulated genes - geometric mean of the 
  #down-regulated genes in the signature
  
  #pull out positive (up-regulated) and negative (down_regulated) 
  #genes from signature
  pos.genes = filter_object$posGeneNames
  neg.genes = filter_object$negGeneNames
  
  #intersect these with the genes present in the expression matrix
  genes_present = row.names(datasetObject$expr)
  pos.genes_present = intersect(pos.genes, genes_present)
  neg.genes_present = intersect(neg.genes, genes_present)
  
  posScore = compute_sample_scores(datasetObject, pos.genes_present)
  negScore = compute_sample_scores(datasetObject, neg.genes_present)
  
  totalScore = posScore - negScore
  
  if (sum(abs(totalScore)) != 0) totalScore = as.numeric(scale(totalScore))
  
  if (!suppressMessages) {
    cat("Used ", length(pos.genes_present), 
        "of ", length(pos.genes), " pos genes, and ", length(neg.genes_present), " of ", length(neg.genes), 
        " neg genes \n")
  }
  return(totalScore)
}


compute_sample_scores = function(datasetObject, genes){
  #compute the geometric mean of a group of genes in a 
  #datasetObject consistent with the meta-integrator package
  genes_idx = which(row.names(datasetObject$expr) %in% genes)
  if (length(genes_idx) == 0)   return(rep(0, ncol(datasetObject$expr)))
  
  gene_expr = datasetObject$expr[genes_idx,] 
  #in case only one gene is present, take the vector itself
  if (is.null(nrow(gene_expr))) sample_scores = gene_expr
  if (!is.null(nrow(gene_expr))) sample_scores = apply(gene_expr, 2, geom_mean)
  
  return(sample_scores)
}

geom_mean = function(x){
  
  return(exp(mean(log(x))))
  
}


auroc = function(bool, score) {
  if (length(unique(score)) == 1){
    print('only one value of score: setting AUC to NA')
    return(NA)
  } 
  n1 = sum(!bool)
  n2 = sum(bool)
  U  = sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}


compute_multi_objective_vector = function(filter_object, 
                                          split_dev = F){
  
  #compute vector of objectives that will be optimized
  #all the objectives are defined in terms of AUC's resulting
  #from the score obtained using the delta of geometric means in the signatures (as in the
  #meta-integrator package)
  #the objectives are defined to be in the range [0, 1], where 1 is the best
  #and 0 is the worse 
  
  COVID_scores = lapply(COVID_contrasts, function(x) data.frame(score = calculate_score_fast(filter_object, x), 
                                                                group = x$class))
  COVID_aucs = sapply(COVID_scores, function(x) auroc(x$group, x$score))
  
  COVID_objective_vector = tapply(COVID_aucs, COVID_contrasts_dict$superclass, 
                                  function(x) min(x, na.rm = T))
  
  #in split_dev = TRUE, the multiobjective vector is computed
  #separately for discovery and development set
  if (split_dev){
    
    COVID_objective_vector = tapply(COVID_aucs, interaction(COVID_contrasts_dict$superclass, COVID_contrasts_dict$use), 
                                    function(x) min(x, na.rm = T))
  }
  
  
  non_COVID_scores = lapply(non_COVID_contrasts, function(x) data.frame(score = calculate_score_fast(filter_object, x), 
                                                                        group = x$class))
  non_COVID_aucs = sapply(non_COVID_scores, function(x) auroc(x$group, x$score))
  
  #the definition of objective below penalizes AUC's larger than 0.5 in non COVID contrasts
  #the factor 2 is meant to adjust the range in [0, 1]
  
  non_COVID_objective_vector = tapply(non_COVID_aucs, non_COVID_contrasts_dict$superclass,
                                      function(x) 1 - 2*max(ifelse(x > 0.5, x - 0.5, 0), na.rm = T))
  
  #in split_dev = TRUE, the multiobjective vector is computed
  #separately for discovery and development set
  if (split_dev){
    
    non_COVID_objective_vector = tapply(non_COVID_aucs, interaction(non_COVID_contrasts_dict$superclass, non_COVID_contrasts_dict$use),
                                        function(x) 1 - 2*max(ifelse(x > 0.5, x - 0.5, 0)))
    
  }
  
  
  return(c(COVID_objective_vector, non_COVID_objective_vector))
  
}


scalarize_multi_objective_vector = function(multi_objective_vector, 
                                            weight_COVID_vs_healthy = 1, 
                                            weight_COVID_vs_infection = 1, 
                                            weight_Non.Microbial = 1,
                                            weight_Other = 1, 
                                            weight_Resp = 1){
  #given a multi-objective vector and a set of weights,
  #it scalarizes the multi-objective problem either 
  #through linear weighting 
  weight_vector = c(weight_COVID_vs_healthy, weight_COVID_vs_infection, 
                    weight_Non.Microbial, weight_Other, weight_Resp)
  
  #explicit calculation of linear weighting
  scalar_objective = weight_COVID_vs_healthy * multi_objective_vector['COVID_vs_healhty'] + 
    weight_COVID_vs_infection * multi_objective_vector['COVID_vs_infection'] + 
    weight_Non.Microbial * multi_objective_vector['non.infectious'] + 
    weight_Other * multi_objective_vector['Other'] + 
    weight_Resp * multi_objective_vector['Resp']
  
  #rescale so that the scalarized objective is also within the range [0, 1]
  scalar_objective = as.numeric(scalar_objective/sum(weight_vector))
  
  return(scalar_objective)
  
}

check_signature = function(n_positive, n_negative, max_signature_size){
  
  #penalty for imbalance
  imbalance = max(n_positive, n_negative)/min(n_positive, n_negative)
  if (imbalance > 1.5){
    #print('imbalanced signature, rejected')
    check = FALSE
    return(check)
  } 
  
  #penalty for signature size
  signature_size = n_positive + n_negative
  if (signature_size > max_signature_size){
    #print('signature size > 25, rejected')
    check = FALSE
    return(check)
  } 
  
  check = TRUE
  return(check)
  
}


get_multi_objective_vector_from_GA_object = function(GA_object, filter_object,
                                                     COVID_contrasts, 
                                                     COVID_contrasts_dict, 
                                                     non_COVID_contrasts, 
                                                     non_COVID_contrasts_dict){
  #pool out from a GA solution the corresponding vectors of multi-objectives
  #both from the discovery and from the development studies
  #the result of this function will be used to compute the Pareto front
  genes_in_signature = sort(pool_of_genes[which(GA_object@solution[1, ] == 1)])
  
  #now evaluate signature on all studies (discovery and validation)
  #signature = filter_object
  filter_object$posGeneNames = intersect(genes_in_signature, filter_object$posGeneNames)
  filter_object$negGeneNames = intersect(genes_in_signature, filter_object$negGeneNames)
  
  multi_objective_vector = compute_multi_objective_vector(filter_object = filter_object, 
                                                          COVID_contrasts, 
                                                          COVID_contrasts_dict, 
                                                          non_COVID_contrasts, 
                                                          non_COVID_contrasts_dict, 
                                                          split_dev = T)
  
  multi_objective_vector = t(data.frame(discovery = multi_objective_vector[grep('discovery', names(multi_objective_vector))],
                                        development = multi_objective_vector[grep('development', names(multi_objective_vector))]))
  
  colnames(multi_objective_vector) = gsub('.discovery', '', colnames(multi_objective_vector))
  
  return(multi_objective_vector)
  
}

