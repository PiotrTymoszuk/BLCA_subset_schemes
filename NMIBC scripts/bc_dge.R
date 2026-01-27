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
  
  nmibc_bc$anova <- 
    ## anntotation with Entrez ID
    map2(nmibc_bc$anova, 
         nmibc_bc$annotation, 
         left_join, by = "variable") %>% 
    ## significant effects
    map(mutate, 
        regulation = ifelse(p_adjusted < 0.05 & etasq >= 0.06, 
                            "regulated", "ns"), 
        regulation = factor(regulation, c("regulated", "ns")))
  
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
    ## annotation with Entrez ID
    map2(nmibc_bc$test, 
         nmibc_bc$annotation, 
         left_join, 
         by = "variable") %>% 
    ## significant effects in ANOVA
    map2(nmibc_bc$anova_significant, 
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
  
  nmibc_bc$significant <- nmibc_bc$test %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    map(function(x) if(nrow(x) == 0) NULL else x) %>% 
    compact %>% 
    map(blast, clust_id) %>% 
    transpose %>% 
    map(map, blast, regulation) %>% 
    map(transpose) %>% 
    map(map, map, ~.$variable)
  
  ## shared by at least three cohorts
  
  nmibc_bc$cmm_significant <- nmibc_bc$significant %>% 
    map(map, shared_features, m = 3) %>% 
    map(map, as.character)
  
# Number of significant effects ------
  
  insert_msg("Numbers of significant effects")
  
  ## in single cohorts
  
  nmibc_bc$numbers$cohorts <- nmibc_bc$test %>% 
    map(count, clust_id, regulation) %>% 
    map(group_by, clust_id) %>% 
    map(mutate, 
        n_total = sum(n), 
        percent = n/n_total * 100) %>% 
    map(ungroup) %>% 
    map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
    compress(names_to = "cohort")
  
  ## shared significant effects
  
  nmibc_bc$numbers$common <- nmibc_bc$cmm_significant %>% 
    map(map_dbl, length) %>% 
    map(compress, names_to = "regulation", values_to = "n") %>% 
    compress(names_to = "clust_id") %>% 
    mutate(cohort = "common", 
           clust_id = factor(clust_id, globals$cluster_levels), 
           regulation = factor(regulation, globals$regulation_levels), 
           n_total = length(reduce(nmibc_bc$variables, union)), 
           percent = n/n_total * 100)
  
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