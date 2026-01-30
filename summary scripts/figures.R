# Figures with the results of comparison of transcriptomes of the 
# bladder cancer clusters, UROMOL classes of NMIBC, and consensus classes of 
# MIBC.
#
# Functional code generated in this script will be used in the main bladder 
# cancer cluster repository.

  insert_head()
  
# container -----
  
  sum_figs <- list()
  
# figures ---------
  
  insert_msg("Figures")
  
  ## upper panel: similarity and shared markers: 
  ## graph plots
  
  sum_figs$upper <- 
    list(graphs = sum_gplots$graph_plots$j_n_regulated, 
         markers = sum_teplots$plots %>% 
           map(~.x + 
                 expand_limits(x = -0.05) + 
                 expand_limits(x = 1.05))) %>% 
    transpose %>% 
    map2(., 
         list(c("NMIBC subsets similarity", 
                "NMIBC shared markers"), 
              c("MIBC subsets similarity", 
                "MIBC shared markers")), 
         function(ent, title) map2(ent, title, 
                                   ~.x + 
                                     labs(title = .y) + 
                                     theme(legend.position = "none"))) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   rel_widths = c(1.2, 0.8), 
                   align = "hv", 
                   axis = "tblr", 
                   labels = c("A", "B"), 
                   label_size = 10))
  
  ## bottom panels: quantification and function of gene overlaps
  
  sum_figs$bottom <- 
    list(nmibc = 
           sum_veplots$annotated_plots$nmibc[c("bc.#1|uromol.2b.upregulated", 
                                               "bc.#2|uromol.2b.upregulated", 
                                               "bc.#3|uromol.2a.downregulated", 
                                               "bc.#3|uromol.3.downregulated")], 
         mibc = 
           sum_veplots$annotated_plots$mibc[c("bc.#1|consensus.Stroma-rich.upregulated", 
                                              "bc.#2|consensus.Ba/Sq.upregulated", 
                                              "bc.#3|consensus.LumP.upregulated", 
                                              "bc.#3|consensus.LumP.downregulated")]) %>% 
    map(map, ~.x + theme(legend.position = "none")) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   align = "hv", 
                   axis = "tblr", 
                   labels = c("C", "", "", ""), 
                   label_size = 10))
  
  ## the entire figures
  
  sum_figs[c("nmibc", "mibc")] <- 
    map2(sum_figs$upper, 
         sum_figs$bottom, 
         plot_grid, 
         nrow = 2, 
         rel_heights = c(1, 1.2))
  
  sum_figs[c("nmibc", "mibc")] <- sum_figs[c("nmibc", "mibc")] %>% 
    list(label = c("transcriptome_comparison_clusters_uromol_classes_nmibc", 
                   "transcriptome_comparison_clusters_consensus_mibc"), 
         ref_name = names(.), 
         caption = paste("Comparison of differentially regulated transcriptomes", 
                         "of bladder cancer clusters and", 
                         c("UROMOL molecular classes", 
                           "consensus molecular classes"), 
                         c("in non-muscle invasive bladder cancer (NMIBC).", 
                           "in muscle invasive bladder cancer (MIBC)."))) %>% 
    pmap(as_figure, 
         w = 190, 
         h = 210)
  
# Expression of selected markers ---------
  
  insert_msg("Expression of selected markers")

  ## NMIBC: VIM and IL32 in the E-MTAB4321 cohort
  ## MIBC: SMOC2, IL1RAP, SLC44A3 in the TCGA cohort
  
  sum_figs[c("nmibc_markers", "mibc_markers")] <- 
    list(sum_deplots$plots$nmibc$emtab4321 %>% 
           map(~.x[c("VIM", "IL32")]), 
         sum_deplots$plots$mibc$tcga %>% 
           map(~.x[c("SMOC2", "IL1RAP", "SLC44A3")])) %>% 
    map(transpose) %>% 
    map(unlist, recursive = FALSE) %>% 
    map(map, ~.x + theme(legend.position = "none")) %>% 
    map2(., c(1.2, 1.5), 
         ~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   rel_widths = c(1, .y), 
                   align = "hv", 
                   axis = "tblr"))
    
  sum_figs[c("nmibc_markers", "mibc_markers")] <- 
    sum_figs[c("nmibc_markers", "mibc_markers")] %>% 
    list(label = c("selected_markers_bladder_cancer_clusters_uromol_nmibc", 
                   "selected_markers_bladder_cancer_clusters_consensus_mibc"), 
         ref_name = names(.), 
         caption = paste("Expression of selected marker genes shared between", 
                         c("bladder cancer clusters and UROMOL molecular classes", 
                           "bladder cancer clusters and consensus molecular classes"), 
                         c("in NMIBC.", "in MIBC.")), 
         w = c(160, 180), 
         h = c(140, 210)) %>% 
    pmap(as_figure)

# Saving the figures on the disc ---------
  
  insert_msg("Saving the figures on the disc")
  
  sum_figs[c("upper", "bottom")] <- NULL
  
  sum_figs <- compact(sum_figs)
  
  sum_figs %>% 
    walk(pickle, 
         path = "./report/figures", 
         format = "pdf", 
         device = cairo_pdf)
    
# END -----
  
  insert_tail()