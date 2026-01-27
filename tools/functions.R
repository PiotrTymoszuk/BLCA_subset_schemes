# Accessory functions for the project

# differential expression analyses ---------

  format_anova <- function(anova, 
                           annotation, 
                           p_cutoff = 0.05, 
                           etasq_cutoff = 0.06) {
    
    ## consistent formatting of one-way ANOVA results: 
    ## appending with Entrez ID and identificaiton of 
    ## significant effects with cutoffs of adjusted p and effect size
    
    anova %>% 
      left_join(annotation, by = "variable") %>% 
      mutate(regulation = ifelse(p_adjusted < p_cutoff & etasq >= etasq_cutoff, 
                                 "regulated", "ns"), 
             regulation = factor(regulation, c("regulated", "ns")))

  }
  
  format_posthoc <- function(test, 
                             annotation, 
                             anova_significant, 
                             p_cutoff = 0.05, 
                             regulation_levels = globals$regulation_levels) {
    
    ## consistent formatting of post-hoc one-sample T test for 
    ## effects in a subset as compared with the cohort average
    ## annotation with Entrez ID and identification of significant 
    ## effects with a text vector of variables found significant in one-way 
    ## ANOVA and a cutoff of adjusted p value
    
    test %>% 
      left_join(annotation, 
                by = "variable") %>% 
      mutate(anova_significant = ifelse(variable %in% anova_significant, 
                                        "yes", "no"), 
             anova_significant = factor(anova_significant, c("no", "yes"))) %>% 
      mutate(regulation = ifelse(anova_significant == "no" | p_adjusted >= 0.05, 
                                 "ns", 
                                 ifelse(deviation_center > 0, 
                                        "upregulated", 
                                        "downregulated")), 
             regulation = factor(regulation, regulation_levels))
    
  }
  
  find_significant <- function(test_lst, split_fct) {
    
    ## extracts significant effects from post-hoc test results in a list
    
    test_lst %>% 
      map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
      map(function(x) if(nrow(x) == 0) NULL else x) %>% 
      compact %>% 
      map(blast, all_of(split_fct)) %>% 
      transpose %>% 
      map(compact) %>% 
      map(map, blast, regulation) %>% 
      map(transpose) %>% 
      map(map, map, ~.$variable)
    
  }
  
  count_significant <- function(test_lst, split_fct) {
    
    ## counts significant effects in a list of post-hoc results
    
    test_lst %>% 
      map(count, .data[[split_fct]], regulation) %>% 
      map(group_by, .data[[split_fct]]) %>% 
      map(mutate, 
          n_total = sum(n), 
          percent = n/n_total * 100) %>% 
      map(ungroup) %>% 
      map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
      compress(names_to = "cohort")
    
  }
  
  count_cmm_significant <- function(x, 
                                    split_fct, 
                                    universe, 
                                    split_levels, 
                                    regulation_levels = globals$regulation_levels) {
    
    ## counts common significant features
    
    x %>% 
      map(map_dbl, length) %>% 
      map(compress, names_to = "regulation", values_to = "n") %>% 
      compress(names_to = split_fct) %>% 
      mutate(cohort = "common", 
             !!split_fct := factor(.data[[split_fct]], split_levels), 
             regulation = factor(regulation, regulation_levels), 
             n_total = length(universe), 
             percent = n/n_total * 100)
    
  }
  
  find_overlaps <- function(x, 
                            pairs, 
                            regulation_levels = globals$regulation_levels) {
    
    ## extracts shared features from a list of character vectors (x), 
    ## the overlaps to check are defined by pairs: a list with vectors of 
    ## names of elements in x
    
    ## pair names and identifiers
    
    pair_names <- map(pairs, paste, collapse = "|")
    
    pairs <- set_names(pairs, pair_names)
    
    pair_df <- map(pairs, 
                   ~tibble(element1 = .x[[1]], 
                           element2 = .x[[2]]))
    
    pair_df <- compress(pair_df, names_to = "pair_id")
    
    ## overlaps
    
    res <- map(pairs, ~x[.x])
    
    res <- map(res, reduce, intersect)
    
    res <- map(res, 
               ~tibble(feature = .x))
    
    res <- compress(res, names_to = "pair_id")
    
    res <- left_join(res, pair_df, by = "pair_id")
    
    res <- mutate(res, 
                  variable = stri_split_fixed(feature, 
                                              pattern = "|", 
                                              simplify = TRUE)[, 1], 
                  regulation = stri_split_fixed(feature, 
                                                pattern = "|", 
                                                simplify = TRUE)[, 2], 
                  regulation = factor(regulation, regulation_levels))
    
    return(res)
    
  }
  
# General utilities --------
  
  mtx2long <- function(x, 
                       row_var = "variable1", 
                       col_var = "variable2", 
                       value_var = "value") {
    
    ## converts a matrix to long-format data frame: 
    ## variable 1: factors in rows, 
    ## variable 2: factors in columns
    
    x <- as.data.frame(x)
    
    cols_names <- colnames(x)
    
    x <- rownames_to_column(x, row_var)
    
    pivot_longer(x, 
                 cols = all_of(cols_names), 
                 names_to = col_var, 
                 values_to = value_var)
    
  }
  
  
# END --------
