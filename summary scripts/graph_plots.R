# Plots for similarity graphs of bladder cancer clusters, UROMOL classes, 
# and consensus classes of MIBC: 
# 
# 1) Heat map plots with similarity coefficients, numbers of shared regulated 
# genes, and text labels with markers and processes.
#
# 2) Graph plots

  insert_head()
  
# container ---------
  
  sum_gplots <- list()
  
# analysis globals ----------
  
  insert_msg("Analysis globals")
  
  ## lexicons of the molecular subsets
  
  sum_gplots$subset_lexicon <- 
    list(x = globals[c("uromol_levels", "consensus_levels")], 
         y = c("NMIBC_class", "consensusClass"), 
         z = c("uromol", "consensus")) %>% 
    pmap(function(x, y, z) globals$cluster_levels %>% 
           map_dfr(~tibble(clust_id = factor(.x, globals$cluster_levels), 
                           !!y := x, 
                           pair_id = paste0("bc.", clust_id, 
                                            "|", z, ".", 
                                            .data[[y]])))) %>% 
    set_names(c("nmibc", "mibc"))
  
  sum_gplots$subset_lexicon$nmibc <- 
    sum_gplots$subset_lexicon$nmibc %>% 
    mutate(NMIBC_class = factor(NMIBC_class, globals$uromol_levels))
  
  sum_gplots$subset_lexicon$mibc <- 
    sum_gplots$subset_lexicon$mibc %>% 
    mutate(consensusClass = factor(consensusClass, globals$consensus_levels))

  ## graph objects
  
  sum_gplots$graph_obj <- 
    list(nmibc = nmibc_shared, mibc = mibc_shared) %>% 
    map(~.x$simil_graph)
  
  ## data frames for heat maps
  
  sum_gplots$hm_data <- sum_gplots$graph_obj %>% 
    map(edge.attributes) %>% 
    map(as_tibble) %>% 
    ## total regulated genes
    map(mutate, 
        n_regulated = downregulated + upregulated, 
        n_label = paste0("up: n = ", upregulated, 
                         "\ndown: n = ", downregulated))
  
  sum_gplots$hm_data <- 
    map2(sum_gplots$hm_data, 
         sum_gplots$subset_lexicon, 
         left_join, by = "pair_id")
  
  sum_gplots$subset_vars <- 
    set_names(globals$system_levels[2:3], 
              c("nmibc", "mibc"))
  
  sum_gplots$subset_labels <- 
    set_names(globals$system_labels[c("NMIBC_class", "consensusClass")], 
              c("nmibc", "mibc"))
  
  ## common arrow
  
  sum_gplots$cmm_arrow <- arrow(angle = 15, 
                                length = unit(0.5, "lines"), 
                                type = "closed")
  
# Heat maps of edge attributes ------
  
  insert_msg("Heat maps")
  
  ## J coefficients, numbers of regulated, down- and upregulated genes
  
  for(i in names(sum_gplots$hm_data)) {
    
    sum_gplots$attr_hm_plots[[i]] <- 
      list(x = c("j", "n_regulated", "downregulated", "upregulated"), 
           y = c("Similarity", 
                 "Common regulated genes", 
                 "Common downregulated genes", 
                 "Common upregulated genes"), 
           z = c("J", rep("# genes", 3))) %>% 
      pmap(function(x, y, z) sum_gplots$hm_data[[i]] %>% 
             ggplot(aes(x = clust_id, 
                        y = .data[[sum_gplots$subset_vars[[i]]]], 
                        fill = .data[[x]])) + 
             geom_tile(color = "black") + 
             scale_fill_gradient2(low = globals$regulation_colors[["downregulated"]], 
                                  high = globals$regulation_colors[["upregulated"]], 
                                  mid = "white", 
                                  midpoint = sum_gplots$hm_data[[i]][[x]] %>% 
                                    range(na.rm = TRUE) %>% 
                                    mean) + 
             globals$common_theme + 
             labs(title = paste(y, toupper(i), sep = ", "), 
                  fill = z, 
                  x = globals$system_labels[["clust_id"]], 
                  y = sum_gplots$subset_labels[[i]])) %>% 
      set_names(c("j", "n_regulated", "downregulated", "upregulated"))
    
    ## heat map plots with fill color corresponding to J and 
    ## text labels with 
    ## numbers of up- and downregulated genes, the top markers 
    ## and biological processes
    
    sum_gplots$attr_hm_plots[[i]][["markers"]] <- 
      sum_gplots$attr_hm_plots[[i]][["j"]] + 
      geom_text(aes(label = ifelse(j >= 0.1, marker_label, NA)), 
                size = 2, 
                fontface = "italic") + 
      labs(title = paste("Shared markers,", toupper(i)))
    
    sum_gplots$attr_hm_plots[[i]][["go_downregulated"]] <- 
      sum_gplots$attr_hm_plots[[i]][["j"]] + 
      geom_text(aes(label = ifelse(j >= 0.1, label_downregulated, NA)), 
                size = 2) + 
      labs(title = paste("GO enrichment, downregulated genes,", toupper(i)))
    
    sum_gplots$attr_hm_plots[[i]][["go_upregulated"]] <- 
      sum_gplots$attr_hm_plots[[i]][["j"]] + 
      geom_text(aes(label = ifelse(j >= 0.1, label_upregulated, NA)), 
                size = 2) + 
      labs(title = paste("GO enrichment, upregulated genes,", toupper(i)))
    
    sum_gplots$attr_hm_plots[[i]][["j_n_regulated"]] <- 
      sum_gplots$attr_hm_plots[[i]][["j"]] + 
      geom_text(aes(label = n_label), 
                size = 2.5)

    ## text labels for the tiles of heat maps with J 
    ## coefficients and numbers of regulated genes
    
    sum_gplots$attr_hm_plots[[i]][["j"]] <- 
      sum_gplots$attr_hm_plots[[i]][["j"]] + 
      geom_text(aes(label = paste("J =", signif(j, 2))), 
                size = 2.5)
    
    sum_gplots$attr_hm_plots[[i]][["n_regulated"]] <- 
      sum_gplots$attr_hm_plots[[i]][["n_regulated"]] + 
      geom_text(aes(label = n_label), 
                size = 2.5)
    
    sum_gplots$attr_hm_plots[[i]][c("downregulated", 
                                    "upregulated")] <- 
      sum_gplots$attr_hm_plots[[i]][c("downregulated", 
                                      "upregulated")] %>% 
      map2(., names(.), 
           ~.x + 
             geom_text(aes(label = paste("n =", .data[[.y]])), 
                       size = 2.5))
    
  }
  
# Plots of the similarity graphs -----------
  
  insert_msg("Plots of similarity graphs")
  
  ## base plots: vertices and edges, edge color and width codes for 
  ## J
  
  sum_gplots$graph_plots$base <- 
    list(x = sum_gplots$graph_obj, 
         plot_title = toupper(names(sum_gplots$graph_obj))) %>% 
    pmap(plot, 
         layout = layout.fruchterman.reingold, 
         vertex_color = NA, 
         label_vertices = FALSE, 
         weighting_order = 1, 
         cust_theme = globals$net_theme + 
           theme(plot.subtitle = element_blank()), 
         seed = 1324) %>% 
    map(~.x + 
          geom_nodelabel(aes(label = sub_label, 
                             fill = system), 
                         size = 3, 
                         fontface = "bold", 
                         label.padding = unit(0.15, "lines"), 
                         show.legend = FALSE) + 
          scale_fill_manual(values = globals$system_colors, 
                            labels = globals$system_labels, 
                            name = "classification\nsystem") + 
          scale_linewidth(limits = c(0, 0.4), 
                          range = c(0.15, 2), 
                          name = "similarity, J") + 
          scale_alpha_continuous(limits = c(0, 0.4), 
                                 range = c(0.15, 1), 
                                 name = "similarity, J"))

  ## plots with edge labels with J and numbers of regulated genes
  
  sum_gplots$graph_plots$j_n_regulated <- 
    list(x = sum_gplots$graph_plots$base, 
         y = c(0.1, 0), 
         z = c(-0.15, -0.15)) %>% 
    pmap(function(x, y, z) x  + 
           geom_edgelabel_repel(aes(label = ifelse(j >= 0.1, 
                                                   edge_lab, NA)), 
                                color = "plum4", 
                                size = 2.5, 
                                box.padding = unit(2, "lines"),
                                label.padding = unit(0.1, "lines"), 
                                label.size = 0.1, 
                                force = 4, 
                                nudge_x = y, 
                                nudge_y = z, 
                                arrow = sum_gplots$cmm_arrow))
  
  ## plots with edge labels with the shared markers
  
  sum_gplots$graph_plots$markers <- 
    list(x = sum_gplots$graph_plots$base, 
         y = c(0.1, 0), 
         z = c(-0.15, -0.15)) %>% 
    pmap(function(x, y, z) x  + 
           geom_edgelabel_repel(aes(label = ifelse(j >= 0.1, marker_label, NA)), 
                                color = "aquamarine4", 
                                fontface = "italic", 
                                size = 2, 
                                box.padding = unit(4, "lines"),
                                label.padding = unit(0.1, "lines"), 
                                label.size = 0.1, 
                                force = 4, 
                                nudge_x = y, 
                                nudge_y = z, 
                                arrow = sum_gplots$cmm_arrow))
  
  ## plots with edge labels with the shared biological processes
  
  sum_gplots$graph_plots$go_upregulated <- 
    list(x = sum_gplots$graph_plots$base, 
         y = c(0.1, 0), 
         z = c(-0.15, -0.15)) %>% 
    pmap(function(x, y, z) x  + 
           geom_edgelabel_repel(aes(label = ifelse(j >= 0.1, 
                                                   label_upregulated, NA)), 
                                color = "firebrick4", 
                                size = 2.3, 
                                box.padding = unit(4, "lines"),
                                label.padding = unit(0.1, "lines"), 
                                label.size = 0.1, 
                                force = 4, 
                                nudge_x = y, 
                                nudge_y = z, 
                                arrow = sum_gplots$cmm_arrow))
  
  sum_gplots$graph_plots$go_downregulated <- 
    list(x = sum_gplots$graph_plots$base, 
         y = c(0.1, 0), 
         z = c(-0.15, -0.15)) %>% 
    pmap(function(x, y, z) x  + 
           geom_edgelabel_repel(aes(label = ifelse(j >= 0.1, 
                                                   label_downregulated, NA)), 
                                color = "steelblue4", 
                                size = 2.3, 
                                box.padding = unit(4, "lines"),
                                label.padding = unit(0.1, "lines"), 
                                label.size = 0.1, 
                                force = 4, 
                                nudge_x = y, 
                                nudge_y = z, 
                                arrow = sum_gplots$cmm_arrow))
  
# END --------
  
  sum_gplots <- sum_gplots[c("cmm_arrow", "attr_hm_plots", "graph_plots")]
  
  insert_tail()