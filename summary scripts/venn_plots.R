# Ven plots with numbers of concordantly and discordantly regulated genes 
# in bladder cancer clusters, UROMOL classes of NMIBC, and consensus classes 
# of MIBC

  insert_head()
  
# container --------
  
  sum_veplots <- list()
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## numbers of genes 
  
  sum_veplots$numbers <- 
    list(nmibc = nmibc_shared, mibc = mibc_shared) %>% 
    map(~.x$numbers)
  
  sum_veplots$numbers$mibc <- sum_veplots$numbers$mibc %>% 
    mutate(consensusClass = fct_recode(consensusClass, 
                                       `Stroma\nrich` = "Stroma-rich"))
  
  ## verbose labels for the intersecting genes derived from 
  ## GO enrichment analyses
  
  sum_veplots$go_labels <- 
    list(nmibc = nmibc_shared, mibc = mibc_shared) %>% 
    map(~.x$shared_labels)
  
# plots --------
  
  insert_msg("Plots")
  
  sum_veplots$plots <- 
    list(x = sum_veplots$numbers, 
         set2_n_var = c("n_NMIBC_class", 
                        "n_consensusClass"), 
         set2_label = c("NMIBC_class", 
                        "consensusClass")) %>% 
    pmap(euler_from_df, 
         set1_n_var = "n_clust_id", 
         inter_n_var = "n", 
         set1_label = "clust_id", 
         suffix_label = "regulation", 
         plot_names = "pair_id") %>% 
    map(map,
        ~.x + 
          labs(title = .x$labels$title %>% 
                 stri_replace(fixed = "\n", replacement = "-") %>% 
                 paste0(", genes")))
  
# appending the plots with text labels --------
  
  insert_msg("Appending the plots with text labels")
  
  for(i in names(sum_veplots$plots)) {
    
    sum_veplots$annotated_plots[[i]] <- 
      list(x = sum_veplots$plots[[i]][names(sum_veplots$go_labels[[i]])], 
           label = paste("shared processes:", 
                         sum_veplots$go_labels[[i]], 
                         sep = "\n\n")) %>% 
      pmap(euler_add_label, 
           x_offset = 6, 
           txt_size = 2.5) %>% 
      map(~.x + expand_limits(x = 72))
    
  }

# END ---------
  
  sum_veplots <- sum_veplots[c("plots", "annotated_plots")]
  
  insert_tail()