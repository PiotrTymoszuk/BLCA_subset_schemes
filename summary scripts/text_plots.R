# Visualizations of text labels: at the moment only marker genes shared between 
# bladder cancer clusters, UROMOL classes, and consensus classes

  insert_head()
  
# container --------
  
  sum_teplots <- list()
  
# the plotting data --------
  
  insert_msg("The plotting data")
  
  ## data frames with pair IDs and molecular subsets
  ## used later in formatting of GO enrichment and gene counting results
  
  sum_teplots$pair_df <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$simil_test) %>% 
    map(select, 
        pair_id, 
        clust_id, 
        any_of(c("NMIBC_class", "consensusClass"))) %>% 
    map(filter, !duplicated(pair_id))
    
  ## the marker labels extracted from the similarity graphs
  
  sum_teplots$data <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$simil_graph) %>% 
    map(edge.attributes) %>% 
    map(as_tibble) %>% 
    map(select, pair_id, marker_label) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map2(sum_teplots$pair_df, 
         left_join, by = "pair_id") %>% 
    map(relocate, 
        clust_id, 
        any_of(c("NMIBC_class", "consensusClass")))
  
  for(i in names(sum_teplots$data)) {
    
    sum_teplots$data[[i]][1:2] <- sum_teplots$data[[i]][1:2] %>% 
      map_dfc(fct_rev)
    
  }
  
  ## graph representation
  
  sum_teplots$attr_df <- globals$attr_df %>% 
    mutate(type = system == "clust_id")
  
  sum_teplots$graph_obj <- sum_teplots$data %>% 
    map(graph_from_data_frame, directed = FALSE) %>% 
    map(set_vertex_attributes, 
        sum_teplots$attr_df)
  
# grph plots ---------
  
  insert_msg("Graph plots")
  
  ## base plots
  
  sum_teplots$plots <- 
    list(x = sum_teplots$graph_obj, 
         plot_title = paste("Shared markers,", 
                   toupper(names(sum_teplots$data)))) %>% 
    pmap(plot, 
         layout = layout.bipartite, 
         vertex_color = NA, 
         edge_color = NA, 
         cust_theme = globals$net_theme + 
           theme(plot.subtitle = element_blank(), 
                 legend.position = "none")) %>% 
    map(~.x + 
          scale_x_continuous(trans = "reverse") + 
          scale_y_continuous(limits = c(-0.1, 1.1)) + 
          coord_flip())
  
  ## adding the vertices and edges
  
  sum_teplots$plots <- sum_teplots$plots %>% 
    map(~.x + 
          geom_edges(color = "black") + 
          geom_nodelabel(aes(label = sub_label, 
                             fill = system), 
                         fontface = "bold", 
                         label.padding = unit(0.15, "lines"), 
                         size = 2.75) + 
          geom_edgelabel(aes(label = marker_label), 
                         color = "black", 
                         size = 2.1, 
                         fontface = "italic") + 
          scale_fill_manual(values = globals$system_colors, 
                            name = "classification system"))
  
# END --------
  
  sum_teplots <- sum_teplots[c("plots")]
  
  insert_tail()
  