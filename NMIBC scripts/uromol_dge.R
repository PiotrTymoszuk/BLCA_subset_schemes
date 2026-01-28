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
  
  nmibc_uro$anova <- nmibc_uro[c("anova", "annotation")] %>% 
    pmap(format_anova)
  
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
    nmibc_uro[c("test", "annotation", "anova_significant")] %>% 
    pmap(format_posthoc)

# Significant effects ---------
  
  insert_msg("Significant effects")
  
  ## in single cohorts
  
  nmibc_uro$significant <- nmibc_uro$test %>% 
    find_significant(split_fct = "NMIBC_class")
  
  ## shared by at least three cohorts
  
  nmibc_uro$cmm_significant <- nmibc_uro$significant %>% 
    map(map, shared_features, m = 3) %>% 
    map(map, as.character)
  
# Number of significant effects ------
  
  insert_msg("Numbers of significant effects")
  
  ## in single cohorts
  
  nmibc_uro$numbers$cohorts <- nmibc_uro$test %>% 
    count_significant(split_fct = "NMIBC_class")
  
  ## shared significant effects
  
  nmibc_uro$numbers$common <- nmibc_uro$cmm_significant %>% 
    count_cmm_significant(split_fct = "NMIBC_class", 
                          universe = reduce(nmibc_uro$variables, union), 
                          split_levels = globals$uromol_levels)

  nmibc_uro$numbers <- nmibc_uro$numbers %>% 
    reduce(rbind)
  
# consistency of the results: overlaps of differentially regulated genes -------
  
  insert_msg("Consistency of the results")
  
  nmibc_uro$cohort_simil <- nmibc_uro$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(~.x$variable) %>% 
    set_similarity(method = "jaccard")
  
# Markers of the clusters in single cohorts and shared ones --------
  
  insert_msg("Markers of the clusters")
  
  ## the markers selected among ANOVA-significant features 
  ## with AUC >= 0.714
  
  nmibc_uro$roc <- nmibc_uro[c("data", "anova_significant")] %>% 
    set_names(c("data", "variables")) %>% 
    pmap(classify, 
         split_fct = "NMIBC_class") %>% 
    map(~.x$classification)
  
  nmibc_uro$roc <- nmibc_uro[c("roc", "annotation")] %>% 
    pmap(format_roc)
  
  ## markers in single cohorts and shared by at least three cohorts
  
  nmibc_uro$markers <- nmibc_uro$roc %>% 
    find_markers(split_fct = "NMIBC_class")
  
  nmibc_uro$cmm_markers <- nmibc_uro$markers %>% 
    map(shared_features, m = 3) %>% 
    map(as.character)
  
  ## top 100 markers shared by at least three cohorts
  
  nmibc_uro$cmm_top_markers <- nmibc_uro$markers %>% 
    map(map, ~.x[1:100]) %>% 
    map(shared_features, m = 3) %>% 
    map(as.character)
  
# Caching -------
  
  insert_msg("Caching")
  
  nmibc_uro <- 
    nmibc_uro[c("stats", "anova", "test", 
               "significant", "cmm_significant", 
               "numbers", "cohort_simil", 
               "roc", "markers", "cmm_markers", "cmm_top_markers")]
  
  save(nmibc_uro, file = "./cache/nmibc_uro.RData")
  
# END --------
  
  insert_tail()