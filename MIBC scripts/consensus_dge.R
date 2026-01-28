# Differential gene expression for consensus classes of MIBC as compared with 
# the cohort average. One-way ANOVA and one-sample post-hoc T test. 
# Significant effects: pFDR ANOVA < 0.05, eta-square >= 0.06, pFDR T test < 0.05.
#
# Common differentially regulated genes: shared by at least three cohorts.

  insert_head()
  
# container --------
  
  mibc_cons <- list()
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## analysis data: gene expression and consensus class assignment
  
  mibc_cons$assignment <- ex_mibc$assignment %>% 
    map(select, sample_id, consensusClass)
  
  mibc_cons$data <- 
    map2(ex_mibc$assignment, 
         ex_mibc$expression, 
         inner_join, by = "sample_id")
  
  ## analysis data: in the GSE128192 cohort, there's only one NE-like tumor
  ## and hence removed to make ANOVA work
  
  mibc_cons$data$gse128192 <- mibc_cons$data$gse128192 %>% 
    filter(consensusClass != "NE-like") %>% 
    mutate(consensusClass = droplevels(consensusClass))
  
  ## analysis variables 
  
  mibc_cons$annotation <- ex_mibc$annotation %>% 
    map(mutate, variable = gene_symbol)
  
  mibc_cons$variables <- mibc_cons$annotation %>% 
    map(~.x$gene_symbol)
  
# Descriptive stats --------
  
  insert_msg("Descriptive stats")
  
  mibc_cons$stats <- mibc_cons[c("data", "variables")] %>% 
    pmap(fast_num_stats, 
         split_fct = "consensusClass")
  
# One-way ANOVA ---------
  
  insert_msg("One-ways ANOVA")
  
  mibc_cons$anova <- 
    map2(mibc_cons$data, 
         mibc_cons$variables, 
         ~f_one_anova(.x[.y], 
                      f = .x[["consensusClass"]], 
                      as_data_frame = TRUE, 
                      adj_method = "BH")) %>% 
    map(as_tibble)
  
  ## formatting the test results, significant differences 
  ## between the classes
  
  mibc_cons$anova <- mibc_cons[c("anova", "annotation")] %>% 
    pmap(format_anova)
  
  mibc_cons$anova_significant <- mibc_cons$anova %>% 
    map(filter, regulation == "regulated") %>% 
    map(~.x$variable)
  
# Post-hoc T test ----------
  
  insert_msg("Post-hoc T test")

  mibc_cons$test <- mibc_cons[c("data", "variables")] %>% 
    pmap(avg_deviation, 
         split_fct = "consensusClass")
  
  ## formatting the testing results
  
  mibc_cons$test <- 
    mibc_cons[c("test", "annotation", "anova_significant")] %>% 
    pmap(format_posthoc)
  
# Significant effects ---------
  
  insert_msg("Significant effects")
  
  ## in single cohorts
  
  mibc_cons$significant <- mibc_cons$test %>% 
    find_significant(split_fct = "consensusClass")
  
  ## shared by at least three cohorts
  
  mibc_cons$cmm_significant <- mibc_cons$significant %>% 
    map(map, shared_features, m = 3) %>% 
    map(map, as.character)
  
# Number of significant effects ------
  
  insert_msg("Numbers of significant effects")
  
  ## in single cohorts
  
  mibc_cons$numbers$cohorts <- mibc_cons$test %>% 
    count_significant(split_fct = "consensusClass")
  
  ## shared significant effects
  
  mibc_cons$numbers$common <- mibc_cons$cmm_significant %>% 
    count_cmm_significant(split_fct = "consensusClass", 
                          universe = reduce(mibc_cons$variables, union), 
                          split_levels = globals$consensus_levels)
  
  mibc_cons$numbers <- mibc_cons$numbers %>% 
    reduce(rbind)
  
# consistency of the results: overlaps of differentially regulated genes -------
  
  insert_msg("Consistency of the results")
  
  mibc_cons$cohort_simil <- mibc_cons$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(~.x$variable) %>% 
    set_similarity(method = "jaccard")
  
# Markers of the clusters in single cohorts and shared ones --------
  
  insert_msg("Markers of the clusters")
  
  ## the markers selected among ANOVA-significant features 
  ## with AUC >= 0.714
  
  mibc_cons$roc <- mibc_cons[c("data", "anova_significant")] %>% 
    set_names(c("data", "variables")) %>% 
    pmap(classify, 
         split_fct = "consensusClass") %>% 
    map(~.x$classification)
  
  mibc_cons$roc <- mibc_cons[c("roc", "annotation")] %>% 
    pmap(format_roc)
  
  ## markers in single cohorts and shared by at least three cohorts
  
  mibc_cons$markers <- mibc_cons$roc %>% 
    find_markers(split_fct = "consensusClass")
  
  mibc_cons$cmm_markers <- mibc_cons$markers %>% 
    map(shared_features, m = 3) %>% 
    map(as.character)
  
  ## top 100 markers shared by at least three cohorts
  
  mibc_cons$cmm_top_markers <- mibc_cons$markers %>% 
    map(map, ~.x[1:100]) %>% 
    map(shared_features, m = 3) %>% 
    map(as.character)
  
# Caching -------
  
  insert_msg("Caching")
  
  mibc_cons <- 
    mibc_cons[c("stats", "anova", "test", 
               "significant", "cmm_significant", 
               "numbers", "cohort_simil", 
               "roc", "markers", "cmm_markers", "cmm_top_markers")]
  
  save(mibc_cons, file = "./cache/mibc_cons.RData")
  
# END --------
  
  insert_tail()