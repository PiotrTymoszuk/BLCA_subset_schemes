# Differential gene expression for bladder cancer clusters as compared with 
# the cohort average. One-way ANOVA and one-sample post-hoc T test. 
# Significant effects: pFDR ANOVA < 0.05, eta-square >= 0.06, pFDR T test < 0.05.
#
# Common differentially regulated genes: shared by at least three cohorts.

  insert_head()
  
# container --------
  
  nmibc_bc <- list()
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## analysis data: gene expression and cluster assignment
  
  nmibc_bc$assignment <- ex_nmibc$assignment %>% 
    map(select, sample_id, clust_id)
  
  nmibc_bc$data <- 
    map2(ex_nmibc$assignment, 
         ex_nmibc$expression, 
         inner_join, by = "sample_id")
  
  ## analysis variables 
  
  nmibc_bc$annotation <- ex_nmibc$annotation %>% 
    map(mutate, variable = gene_symbol)
  
  nmibc_bc$variables <- nmibc_bc$annotation %>% 
    map(~.x$gene_symbol)
  
# Descriptive stats --------
  
  insert_msg("Descriptive stats")
  
  nmibc_bc$stats <- nmibc_bc[c("data", "variables")] %>% 
    pmap(fast_num_stats, 
         split_fct = "clust_id")
  
# One-way ANOVA ---------
  
  insert_msg("One-ways ANOVA")
  
  nmibc_bc$anova <- 
    map2(nmibc_bc$data, 
         nmibc_bc$variables, 
         ~f_one_anova(.x[.y], 
                      f = .x[["clust_id"]], 
                      as_data_frame = TRUE, 
                      adj_method = "BH")) %>% 
    map(as_tibble)
  
  ## formatting the test results, significant differences 
  ## between the clusters
  
  nmibc_bc$anova <- nmibc_bc[c("anova", "annotation")] %>% 
    pmap(format_anova)
  
  nmibc_bc$anova_significant <- nmibc_bc$anova %>% 
    map(filter, regulation == "regulated") %>% 
    map(~.x$variable)
  
# Post-hoc T test ----------
  
  insert_msg("Post-hoc T test")

  nmibc_bc$test <- nmibc_bc[c("data", "variables")] %>% 
    pmap(avg_deviation, 
         split_fct = "clust_id")
  
  ## formatting the testing results
  
  nmibc_bc$test <- 
    nmibc_bc[c("test", "annotation", "anova_significant")] %>% 
    pmap(format_posthoc)
  
# Significant effects ---------
  
  insert_msg("Significant effects")
  
  ## in single cohorts
  
  nmibc_bc$significant <- nmibc_bc$test %>% 
    find_significant(split_fct = "clust_id")
  
  ## shared by at least three cohorts
  
  nmibc_bc$cmm_significant <- nmibc_bc$significant %>% 
    map(map, shared_features, m = 3) %>% 
    map(map, as.character)
  
# Number of significant effects ------
  
  insert_msg("Numbers of significant effects")
  
  ## in single cohorts
  
  nmibc_bc$numbers$cohorts <- nmibc_bc$test %>% 
    count_significant(split_fct = "clust_id")
  
  ## shared significant effects
  
  nmibc_bc$numbers$common <- nmibc_bc$cmm_significant %>% 
    count_cmm_significant(split_fct = "clust_id", 
                          universe = reduce(nmibc_bc$variables, union), 
                          split_levels = globals$cluster_levels)
  
  nmibc_bc$numbers <- nmibc_bc$numbers %>% 
    reduce(rbind)
  
# consistency of the results: overlaps of differentially regulated genes -------
  
  insert_msg("Consistency of the results")
  
  nmibc_bc$cohort_simil <- nmibc_bc$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(~.x$variable) %>% 
    set_similarity(method = "jaccard")
  
# Caching -------
  
  insert_msg("Caching")
  
  nmibc_bc <- 
    nmibc_bc[c("stats", "anova", "test", 
               "significant", "cmm_significant", 
               "numbers", "cohort_simil")]
  
  save(nmibc_bc, file = "./cache/nmibc_bc.RData")
  
# END --------
  
  insert_tail()