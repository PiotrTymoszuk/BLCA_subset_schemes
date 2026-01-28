# Plots for the results of differential gene expression analyses: 
# heat map representations for the overlapping differentially regulated 
## genes, and markers of bladder cancer clusters, UROMOL classes and 
# consensus classes of NMIBC

  insert_head()
  
# container -------
  
  sum_deplots <- list()
  
# parallel backend --------
  
  insert_msg("Parallel backend")
  
  plan("multisession")
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## shared differentially expressed genes
  ## shared markers of the subsets
  
  sum_deplots$genes <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$shared_features$gene_symbol) %>% 
    map(unique)
  
  sum_deplots$markers <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$shared_markers) %>% 
    map(unlist) %>% 
    map(unique)
  
  ## assignment and splitting factors
  
  sum_deplots$assignment <- list(nmibc = ex_nmibc, 
                                 mibc = ex_mibc) %>% 
    map(~.x$assignment)
  
  sum_deplots$split_fct <- 
    list(nmibc = c("clust_id", "NMIBC_class"), 
         mibc = c("clust_id", "consensusClass"))
  
  sum_deplots$split_labels <- 
    list(nmibc = globals$system_labels[c("clust_id", "NMIBC_class")], 
         mibc = globals$system_labels[c("clust_id", "consensusClass")])
  
  ## gene expression data
  
  sum_deplots$expression <- 
    list(nmibc = ex_nmibc, mibc = ex_mibc) %>% 
    map(~.x$expression)
  
  for(i in names(sum_deplots$expression)) {
    
    sum_deplots$expression[[i]] <- sum_deplots$expression[[i]] %>% 
      map(select, 
          sample_id, 
          any_of(union(sum_deplots$genes[[i]], 
                       sum_deplots$markers[[i]]))) %>% 
      map2(sum_deplots$assignment[[i]], ., 
           inner_join, by = "sample_id")
    
  }
  
# Heat maps for the overlapping genes, single cohorts --------
  
  insert_msg("Heat maps for the overlapping genes, single cohorts")
  
  for(i in names(sum_deplots$expression)) {
    
    for(j in sum_deplots$split_fct[[i]]) {
      
      sum_deplots$gene_hm_plots[[i]][[j]] <- 
        list(data = sum_deplots$expression[[i]], 
             variables = sum_deplots$expression[[i]] %>% 
               map(names) %>% 
               map(intersect, sum_deplots$genes[[i]]), 
             plot_title = globals$cohort_labs[names(sum_deplots$expression[[i]])] %>% 
               paste(sum_deplots$split_labels[[i]][[j]], sep = ", ")) %>% 
        future_pmap(heat_map, 
                    split_fct = j, 
                    normalize = TRUE, 
                    hide_x_axis_text = TRUE, 
                    x_lab = "cancer sample", 
                    y_lab = "gene", 
                    cust_theme = globals$common_theme + 
                      theme(strip.text.y = element_blank(), 
                            strip.background.y = element_blank(), 
                            axis.text.y = element_blank()), 
                    midpoint = 0, 
                    limits = c(-3, 3), 
                    oob = scales::squish, 
                    .options = furrr_options(seed = TRUE))

    }
    
  }
  
# Heat maps for the overlapping markers, single cohorts ---------
  
  insert_msg("Heat maps for the overlapping markers, single cohorts")
  
  for(i in names(sum_deplots$expression)) {
    
    for(j in sum_deplots$split_fct[[i]]) {
      
      sum_deplots$marker_hm_plots[[i]][[j]] <- 
        list(data = sum_deplots$expression[[i]], 
             variables = sum_deplots$expression[[i]] %>% 
               map(names) %>% 
               map(intersect, sum_deplots$markers[[i]]), 
             plot_title = globals$cohort_labs[names(sum_deplots$expression[[i]])] %>% 
               paste(sum_deplots$split_labels[[i]][[j]], sep = ", ")) %>% 
        future_pmap(heat_map, 
                    split_fct = j, 
                    normalize = TRUE, 
                    hide_x_axis_text = TRUE, 
                    x_lab = "cancer sample", 
                    cust_theme = globals$common_theme + 
                      theme(strip.text.y = element_blank(), 
                            strip.background.y = element_blank(), 
                            axis.title.y = element_blank(), 
                            axis.text.y = element_text(face = "italic")), 
                    midpoint = 0, 
                    limits = c(-3, 3), 
                    oob = scales::squish, 
                    .options = furrr_options(seed = TRUE))
      
      if(i == "mibc") {
        
        sum_deplots$marker_hm_plots[[i]][[j]] <- 
          sum_deplots$marker_hm_plots[[i]][[j]] %>% 
          map(~.x + 
                guides(y = guide_axis(n.dodge = 3)) + 
                theme(axis.text.y = element_text(size = 6)))
        
      }
      
    }
    
  }
  
# Summary heat maps for the overlapping regulated genes, cohort averages --------
  
  insert_msg("Summary heat maps, overlapping genes")
  
  for(i in names(sum_deplots$expression)) {
    
    sum_deplots$summary_gene_hm_plots[[i]] <- 
      list(split_fct = sum_deplots$split_fct[[i]], 
           plot_title = sum_deplots$split_labels[[i]]) %>% 
      pmap(common_heat_map, 
           data = sum_deplots$expression[[i]] %>% 
             map(pad_missing, 
                 variables = sum_deplots$genes[[i]]) %>% 
             set_names(globals$cohort_labs[names(sum_deplots$expression[[i]])]), 
           variables = sum_deplots$genes[[i]], 
           normalize = TRUE, 
           x_lab = "cohort", 
           y_lab = "gene", 
           cust_theme = globals$common_theme + 
             theme(strip.text.y = element_blank(), 
                   strip.background.y = element_blank(), 
                   axis.text.y = element_blank(), 
                   axis.ticks.y = element_blank()), 
           midpoint = 0, 
           limits = c(-1.5, 1.5), 
           oob = scales::squish) %>% 
      set_names(sum_deplots$split_fct[[i]])

  }
  
# Summary heat maps for the overlapping markers, cohort averages --------
  
  insert_msg("Summary heat maps, overlapping markers")
  
  for(i in names(sum_deplots$expression)) {
    
    sum_deplots$summary_marker_hm_plots[[i]] <- 
      list(split_fct = sum_deplots$split_fct[[i]], 
           plot_title = sum_deplots$split_labels[[i]]) %>% 
      pmap(common_heat_map, 
           data = sum_deplots$expression[[i]] %>% 
             map(pad_missing, 
                 variables = sum_deplots$markers[[i]]) %>% 
             set_names(globals$cohort_labs[names(sum_deplots$expression[[i]])]), 
           variables = sum_deplots$markers[[i]], 
           normalize = TRUE, 
           cust_theme = globals$common_theme + 
             theme(strip.text.y = element_blank(), 
                   strip.background.y = element_blank(), 
                   axis.title.y = element_blank(), 
                   axis.title.x = element_blank(), 
                   axis.text.y = element_text(face = "italic")), 
           midpoint = 0, 
           limits = c(-1.5, 1.5), 
           oob = scales::squish) %>% 
      set_names(sum_deplots$split_fct[[i]])
    
    if(i == "mibc") {
      
      sum_deplots$summary_marker_hm_plots[[i]] <- 
        sum_deplots$summary_marker_hm_plots[[i]] %>% 
        map(~.x + 
              guides(y = guide_axis(n.dodge = 3)) + 
              theme(axis.text.y = element_text(size = 6)))
      
    }
    
  }
  
# END ---------
  
  plan("sequential")
  
  insert_tail()