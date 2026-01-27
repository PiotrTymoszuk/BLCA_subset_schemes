# Differential gene expression for UROMOL classes as compared with 
# the cohort average. One-way ANOVA and one-sample post-hoc T test. 
# Significant effects: pFDR ANOVA < 0.05, eta-square >= 0.06, pFDR T test < 0.05.
#
# Common differentially regulated genes: shared by at least three cohorts.

  insert_head()
  
# container --------
  
  nmibc_uro <- list()
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## analysis data: gene expression and class assignment
  
  nmibc_uro$assignment <- ex_nmibc$assignment %>% 
    map(select, sample_id, NMIBC_class)
  
  nmibc_uro$data <- 
    map2(ex_nmibc$assignment, 
         ex_nmibc$expression, 
         inner_join, by = "sample_id")
  
  ## analysis variables 
  
  nmibc_uro$annotation <- ex_nmibc$annotation %>% 
    map(mutate, variable = gene_symbol)
  
  nmibc_uro$variables <- nmibc_uro$annotation %>% 
    map(~.x$gene_symbol)
  
# Descriptive stats --------
  
  insert_msg("Descriptive stats")
  
  nmibc_uro$stats <- nmibc_uro[c("data", "variables")] %>% 
    pmap(fast_num_stats, 
         split_fct = "NMIBC_class")
  
# One-way ANOVA ---------
  
  insert_msg("One-ways ANOVA")
  
  nmibc_uro$anova <- 
    map2(nmibc_uro$data, 
         nmibc_uro$variables, 
         ~f_one_anova(.x[.y], 
                      f = .x[["NMIBC_class"]], 
                      as_data_frame = TRUE, 
                      adj_method = "BH")) %>% 
    map(as_tibble)
  
  ## formatting the test results, significant differences 
  ## between the classes
  
  nmibc_uro$anova <- 
    ## anntotation with Entrez ID
    map2(nmibc_uro$anova, 
         nmibc_uro$annotation, 
         left_join, by = "variable") %>% 
    ## significant effects
    map(mutate, 
        regulation = ifelse(p_adjusted < 0.05 & etasq >= 0.06, 
                            "regulated", "ns"), 
        regulation = factor(regulation, c("regulated", "ns")))
  
  nmibc_uro$anova_significant <- nmibc_uro$anova %>% 
    map(filter, regulation == "regulated") %>% 
    map(~.x$variable)
  
# Post-hoc T test ----------
  
  insert_msg("Post-hoc T test")

  nmibc_uro$test <- nmibc_uro[c("data", "variables")] %>% 
    pmap(avg_deviation, 
         split_fct = "NMIBC_class")
  
  ## formatting the testing results
  
  nmibc_uro$test <-
    ## annotation with Entrez ID
    map2(nmibc_uro$test, 
         nmibc_uro$annotation, 
         left_join, 
         by = "variable") %>% 
    ## significant effects in ANOVA
    map2(nmibc_uro$anova_significant, 
         ~mutate(.x, 
                 anova_significant = ifelse(variable %in% .y, 
                                            "yes", "no"), 
                 anova_significant = factor(anova_significant, c("no", "yes")))) %>% 
    ## significant differences vs cohort average
    map(mutate, 
        regulation = ifelse(anova_significant == "no" | p_adjusted >= 0.05, 
                           "ns", 
                           ifelse(deviation_center > 0, 
                                  "upregulated", 
                                  "downregulated")), 
        regulation = factor(regulation, globals$regulation_levels))
  
# Significant effects ---------
  
  insert_msg("Significant effects")
  
  ## in single cohorts
  
  nmibc_uro$significant <- nmibc_uro$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(function(x) if(nrow(x) == 0) NULL else x) %>% 
    compact %>% 
    map(blast, NMIBC_class) %>% 
    transpose %>% 
    map(map, blast, regulation) %>% 
    map(transpose) %>% 
    map(map, map, ~.$variable)
  
  ## shared by at least three cohorts
  
  nmibc_uro$cmm_significant <- nmibc_uro$significant %>% 
    map(map, shared_features, m = 3) %>% 
    map(map, as.character)
  
# Number of significant effects ------
  
  insert_msg("Numbers of significant effects")
  
  ## in single cohorts
  
  nmibc_uro$numbers$cohorts <- nmibc_uro$test %>% 
    map(count, NMIBC_class, regulation) %>% 
    map(group_by, NMIBC_class) %>% 
    map(mutate, 
        n_total = sum(n), 
        percent = n/n_total * 100) %>% 
    map(ungroup) %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    compress(names_to = "cohort")
  
  ## shared significant effects
  
  nmibc_uro$numbers$common <- nmibc_uro$cmm_significant %>% 
    map(map_dbl, length) %>% 
    map(compress, names_to = "regulation", values_to = "n") %>% 
    compress(names_to = "NMIBC_class") %>% 
    mutate(cohort = "common", 
           NMIBC_class = factor(NMIBC_class, globals$uromol_levels), 
           regulation = factor(regulation, globals$regulation_levels), 
           n_total = length(reduce(nmibc_uro$variables, union)), 
           percent = n/n_total * 100)
  
  nmibc_uro$numbers <- nmibc_uro$numbers %>% 
    reduce(rbind)
  
# consistency of the results: overlaps of differentially regulated genes -------
  
  insert_msg("Consistency of the results")
  
  nmibc_uro$cohort_simil <- nmibc_uro$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(~.x$variable) %>% 
    set_similarity(method = "jaccard")
  
# Caching -------
  
  insert_msg("Caching")
  
  nmibc_uro <- 
    nmibc_uro[c("stats", "anova", "test", 
               "significant", "cmm_significant", 
               "numbers", "cohort_simil")]
  
  save(nmibc_uro, file = "./cache/nmibc_uro.RData")
  
# END --------
  
  insert_tail()