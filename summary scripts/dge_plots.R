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
  
  sum_deplots$shared_markers <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$shared_markers) %>% 
    map(unlist) %>% 
    map(unique)
  
  sum_deplots$distinct_markers <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$distinct_markers) %>% 
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
  
  sum_deplots$split_colors <- 
    globals[c("cluster_colors", "uromol_colors", "consensus_colors")] %>% 
    set_names(c("clust_id", "NMIBC_class", "consensusClass"))

  ## gene expression data
  
  sum_deplots$expression <- 
    list(nmibc = ex_nmibc, mibc = ex_mibc) %>% 
    map(~.x$expression)
  
  for(i in names(sum_deplots$expression)) {
    
    sum_deplots$expression[[i]] <- sum_deplots$expression[[i]] %>% 
      map(select, 
          sample_id, 
          any_of(reduce(list(sum_deplots$genes[[i]], 
                             sum_deplots$shared_markers[[i]], 
                             sum_deplots$distinct_markers[[i]]), 
                        union))) %>% 
      map2(sum_deplots$assignment[[i]], ., 
           inner_join, by = "sample_id")
    
  }

# Heat maps for the overlapping markers, single cohorts ---------
  
  insert_msg("Heat maps for the overlapping markers, single cohorts")
  
  for(i in names(sum_deplots$expression)) {
    
    for(j in sum_deplots$split_fct[[i]]) {
      
      sum_deplots$shared_marker_hm_plots[[i]][[j]] <- 
        list(data = sum_deplots$expression[[i]], 
             variables = sum_deplots$expression[[i]] %>% 
               map(names) %>% 
               map(intersect, sum_deplots$shared_markers[[i]]), 
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
  
      
      sum_deplots$distinct_marker_hm_plots[[i]][[j]] <- 
        list(data = sum_deplots$expression[[i]], 
             variables = sum_deplots$expression[[i]] %>% 
               map(names) %>% 
               map(intersect, sum_deplots$distinct_markers[[i]]), 
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
        
        sum_deplots$shared_marker_hm_plots[[i]][[j]] <- 
          sum_deplots$shared_marker_hm_plots[[i]][[j]] %>% 
          map(~.x + 
                guides(y = guide_axis(n.dodge = 2)) + 
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
    
    sum_deplots$summary_shared_marker_hm_plots[[i]] <- 
      list(split_fct = sum_deplots$split_fct[[i]], 
           plot_title = sum_deplots$split_labels[[i]]) %>% 
      pmap(common_heat_map, 
           data = sum_deplots$expression[[i]] %>% 
             map(pad_missing, 
                 variables = sum_deplots$shared_markers[[i]]) %>% 
             set_names(globals$cohort_labs[names(sum_deplots$expression[[i]])]), 
           variables = sum_deplots$shared_markers[[i]], 
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
    
    sum_deplots$summary_distinct_marker_hm_plots[[i]] <- 
      list(split_fct = sum_deplots$split_fct[[i]], 
           plot_title = sum_deplots$split_labels[[i]]) %>% 
      pmap(common_heat_map, 
           data = sum_deplots$expression[[i]] %>% 
             map(pad_missing, 
                 variables = sum_deplots$distinct_markers[[i]]) %>% 
             set_names(globals$cohort_labs[names(sum_deplots$expression[[i]])]), 
           variables = sum_deplots$distinct_markers[[i]], 
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
              guides(y = guide_axis(n.dodge = 2)) + 
              theme(axis.text.y = element_text(size = 6)))
      
    }
    
  }
  
# Expression of selected markers in single cohorts and molecular subsets ---------
  
  insert_msg("Box plots for single markers")
  
  for(i in names(sum_deplots$expression)) {
    
    for(j in names(sum_deplots$expression[[i]])) {
      
      temp_vars <- c("VIM", "SLA", "IL32", 
                     "DCN", "SMOC2", "SORBS1", 
                     "CD44", "KRT5", "IL1RAP", 
                     "TBX3", "SLC44A3", "SLC14A1", "RNF128") %>% 
        intersect(names(sum_deplots$expression[[i]][[j]]))
      
      for(k in sum_deplots$split_fct[[i]]) {
        
        sum_deplots$plots[[i]][[j]][[k]] <- 
          list(variable = temp_vars, 
              plot_title = paste(html_italic(temp_vars), 
                                  globals$cohort_labs[j], 
                                  sep = ", ")) %>% 
          pmap(plot_variable, 
                      sum_deplots$expression[[i]][[j]], 
                      split_factor = k, 
                      type = "box", 
                      x_n_labs = TRUE, 
                      x_lab = sum_deplots$split_labels[[i]][[k]], 
                      y_lab = expression("log"[2] * " expression"), 
                      cust_theme = globals$common_theme + 
                        theme(plot.title = element_markdown()), 
               point_hjitter = 0) %>% 
          map(~.x + 
                scale_fill_manual(values = sum_deplots$split_colors[[k]])) %>% 
          set_names(temp_vars)
        
      }
      
    }
    
  }

# END ---------

  rm(temp_vars)
  
  sum_deplots$genes <- NULL
  sum_deplots$markers <- NULL
  sum_deplots$assignment <- NULL
  sum_deplots$expression <- NULL

  sum_deplots <- compact(sum_deplots)
  
  plan("sequential")
  
  insert_tail()