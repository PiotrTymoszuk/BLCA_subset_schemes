# Differential gene expression for bladder cancer clusters as compared with 
# the cohort average. One-way ANOVA and one-sample post-hoc T test. 
# Significant effects: pFDR ANOVA < 0.05, eta-square >= 0.06, pFDR T test < 0.05.
#
# Common differentially regulated genes: shared by at least three cohorts.

  insert_head()
  
# container --------
  
  mibc_bc <- list()
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## analysis data: gene expression and cluster assignment
  
  mibc_bc$assignment <- ex_mibc$assignment %>% 
    map(select, sample_id, clust_id)
  
  mibc_bc$data <- 
    map2(ex_mibc$assignment, 
         ex_mibc$expression, 
         inner_join, by = "sample_id")
  
  ## analysis variables 
  
  mibc_bc$annotation <- ex_mibc$annotation %>% 
    map(mutate, variable = gene_symbol)
  
  mibc_bc$variables <- mibc_bc$annotation %>% 
    map(~.x$gene_symbol)
  
# Descriptive stats --------
  
  insert_msg("Descriptive stats")
  
  mibc_bc$stats <- mibc_bc[c("data", "variables")] %>% 
    pmap(fast_num_stats, 
         split_fct = "clust_id")
  
# One-way ANOVA ---------
  
  insert_msg("One-ways ANOVA")
  
  mibc_bc$anova <- 
    map2(mibc_bc$data, 
         mibc_bc$variables, 
         ~f_one_anova(.x[.y], 
                      f = .x[["clust_id"]], 
                      as_data_frame = TRUE, 
                      adj_method = "BH")) %>% 
    map(as_tibble)
  
  ## formatting the test results, significant differences 
  ## between the clusters
  
  mibc_bc$anova <- mibc_bc[c("anova", "annotation")] %>% 
    pmap(format_anova)
  
  mibc_bc$anova_significant <- mibc_bc$anova %>% 
    map(filter, regulation == "regulated") %>% 
    map(~.x$variable)
  
# Post-hoc T test ----------
  
  insert_msg("Post-hoc T test")

  mibc_bc$test <- mibc_bc[c("data", "variables")] %>% 
    pmap(avg_deviation, 
         split_fct = "clust_id")
  
  ## formatting the testing results
  
  mibc_bc$test <- 
    mibc_bc[c("test", "annotation", "anova_significant")] %>% 
    pmap(format_posthoc)
  
# Significant effects ---------
  
  insert_msg("Significant effects")
  
  ## in single cohorts
  
  mibc_bc$significant <- mibc_bc$test %>% 
    find_significant(split_fct = "clust_id")
  
  ## shared by at least three cohorts
  
  mibc_bc$cmm_significant <- mibc_bc$significant %>% 
    map(map, shared_features, m = 3) %>% 
    map(map, as.character)
  
# Number of significant effects ------
  
  insert_msg("Numbers of significant effects")
  
  ## in single cohorts
  
  mibc_bc$numbers$cohorts <- mibc_bc$test %>% 
    count_significant(split_fct = "clust_id")
  
  ## shared significant effects
  
  mibc_bc$numbers$common <- mibc_bc$cmm_significant %>% 
    count_cmm_significant(split_fct = "clust_id", 
                          universe = reduce(mibc_bc$variables, union), 
                          split_levels = globals$cluster_levels)
  
  mibc_bc$numbers <- mibc_bc$numbers %>% 
    reduce(rbind)
  
# consistency of the results: overlaps of differentially regulated genes -------
  
  insert_msg("Consistency of the results")
  
  mibc_bc$cohort_simil <- mibc_bc$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(~.x$variable) %>% 
    set_similarity(method = "jaccard")
  
# Markers of the clusters in single cohorts and shared ones --------
  
  insert_msg("Markers of the clusters")
  
  ## the markers selected among ANOVA-significant features 
  ## with AUC >= 0.714
  
  mibc_bc$roc <- mibc_bc[c("data", "anova_significant")] %>% 
    set_names(c("data", "variables")) %>% 
    pmap(classify, 
         split_fct = "clust_id") %>% 
    map(~.x$classification)
  
  mibc_bc$roc <- mibc_bc[c("roc", "annotation")] %>% 
    pmap(format_roc)
  
  ## markers in single cohorts and shared by at least three cohorts
  
  mibc_bc$markers <- mibc_bc$roc %>% 
    find_markers(split_fct = "clust_id")
  
  mibc_bc$cmm_markers <- mibc_bc$markers %>% 
    map(shared_features, m = 3) %>% 
    map(as.character)
  
  ## top 100 markers shared by at least three cohorts
  
  mibc_bc$cmm_top_markers <- mibc_bc$markers %>% 
    map(map, ~.x[1:100]) %>% 
    map(shared_features, m = 3) %>% 
    map(as.character)
  
# Caching -------
  
  insert_msg("Caching")
  
  mibc_bc <- 
    mibc_bc[c("stats", "anova", "test", 
              "significant", "cmm_significant", 
              "numbers", "cohort_simil", 
              "roc", "markers", "cmm_markers", "cmm_top_markers")]
  
  save(mibc_bc, file = "./cache/mibc_bc.RData")
  
# END --------
  
  insert_tail()