# Similarity and dissimilarity of the bladder cancer clusters and consensus 
# classes in MIBC: 
#
# 1) Jaccard similarities between the clusters and consensus classes in respect 
# to genes found differentially regulated in at least three cohorts.
#
# 2) Identification of shared upregulated and downregulated genes 
# between the bladder cancer clusters and consensus classes
#
# 3) GO enrichment analyses for the overlapping up- and downregulated genes. 
# Significant enrichment considered for at least 3 genes assigned to the term, 
# pFDR < 0.05, and OR < 0.05.
#
# 4) Construction of a weighted similarity graph with Jaccard similarities.

  insert_head()
  
# container -------
  
  mibc_shared <- list()
  
# parallel backend --------
  
  insert_msg("Parallel backend")
  
  plan("multisession")
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## common differentially regulated genes in bladder cancer clusters 
  ## and common regulated genes in UROMOL classes
  
  mibc_shared$genes[c("bc", "consensus")] <- 
    list(mibc_bc, mibc_cons) %>% 
    map(~.x$cmm_significant)
  
  ## attributes for vertices of similarity graphs
  
  mibc_shared$attr_df <- globals$attr_df
  
  ## gene universe, to be used in GO enrichment analyses
  
  mibc_shared$universe <- mibc_bc$anova %>% 
    map(~.x$entrez_id) %>% 
    reduce(union)
  
# Jaccard distances --------
  
  insert_msg("Jaccard distances")
  
  ## input data: vectors of simplified gene regulation information: 
  ## the upregulate and downregulated features are appended with suffixes
  
  mibc_shared$simil_data <- mibc_shared$genes %>% 
    map(map, ~map2(.x, names(.x), paste, sep = "|")) %>% 
    map(map, reduce, union) %>% 
    unlist(recursive = FALSE)
  
  ## Jaccard distances
  
  mibc_shared$simil_test <- mibc_shared$simil_data %>% 
    set_similarity(method = "jaccard")
  
  mibc_shared$simil_test <- 
    mibc_shared$simil_test[c("bc.#1", "bc.#2", "bc.#3"), 
                           c("consensus.Stroma-rich", 
                             "consensus.LumU", 
                             "consensus.LumNS", 
                             "consensus.Ba/Sq", 
                             "consensus.LumP", 
                             "consensus.NE-like")] %>% 
    mtx2long(row_var = "clust_id", 
             col_var = "consensusClass", 
             value_var = "j") %>% 
    mutate(clust_id = stri_split_fixed(clust_id, 
                                       pattern = ".", 
                                       simplify = TRUE)[, 2], 
           clust_id = factor(clust_id, globals$cluster_levels), 
           consensusClass = stri_split_fixed(consensusClass, 
                                             pattern = ".", 
                                             simplify = TRUE)[, 2], 
           consensusClass = factor(consensusClass, globals$consensus_levels), 
           pair_id = paste(paste0("bc.", clust_id), 
                           paste0("consensus.", consensusClass), 
                           sep = "|"))
  
  ## a data frame with pair IDs and molecular subsets
  ## used later in formatting of GO enrichment and gene counting results

  mibc_shared$pair_df <- 
    mibc_shared$simil_test[, c("pair_id", "clust_id", "consensusClass")] %>% 
    filter(!duplicated(pair_id))

# Identification of shared differentially regulated genes ---------
  
  insert_msg("Identification of shared regulated genes")
  
  ## pairs of molecular subsets
  
  mibc_shared$pairs <- paste0("consensus.", globals$consensus_levels) %>% 
    map(function(x) paste0("bc.", globals$cluster_levels) %>% 
          map(~c(.x, x))) %>% 
    unlist(recursive = FALSE)
  
  ## common regulated features, appended with Entrez IDs

  mibc_shared$shared_features <- 
    find_overlaps(mibc_shared$simil_data, 
                  mibc_shared$pairs) %>% 
    mutate(gene_symbol = variable, 
           entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              keytype = "SYMBOL", 
                              column = "ENTREZID"), 
           clust_id = stri_split_fixed(element1, 
                                       pattern = ".", 
                                       simplify = TRUE)[, 2], 
           clust_id = factor(clust_id, 
                             globals$cluster_levels), 
           consensusClass = stri_split_fixed(element2, 
                                             pattern = ".", 
                                             simplify = TRUE)[, 2], 
           consensusClass = factor(consensusClass, 
                                   globals$uromol_levels))

# Numbers of genes in the intersections ---------
  
  insert_msg("Numbers of genes in the intersections")
  
  mibc_shared$numbers <- mibc_shared$shared_features %>% 
    count(pair_id, regulation) %>% 
    left_join(mibc_shared$pair_df, by = "pair_id")
  
# GO enrichment analyses for the overlaps with at least 10 genes ----------
  
  insert_msg("GO enrichment analyses for the overlaps")
  
  mibc_shared$go_test <- mibc_shared$shared_features %>% 
    blast(pair_id, regulation) %>% 
    map(~.x$entrez_id) %>% 
    map(function(x) if(length(x) < 10) NULL else x) %>% 
    compact %>% 
    future_map(GOana, 
               universe = mibc_shared$universe, 
               ontology = "BP", 
               .options = furrr_options(seed = 1232))
  
  ## formatting the GO enrichment results
  
  mibc_shared$go_test <- mibc_shared$go_test %>% 
    compress(names_to = "condition") %>% 
    mutate(regulation = stri_extract(condition, 
                                     regex = "upregulated|downregulated"), 
           regulation = factor(regulation, globals$regulation_levels), 
           pair_id = stri_replace(condition, 
                                  regex = "\\.(upregulated|downregulated)", 
                                  replacement = "")) %>% 
    left_join(mibc_shared$pair_df, 
              by = "pair_id")
  
  ## identification of significant GO terms
  
  mibc_shared$go_test <- mibc_shared$go_test %>% 
    filter(n_go_dge >= 3) %>% 
    mutate(enrichment = ifelse(p_adjusted >= 0.05 | or < 1.44, 
                               "ns", "enriched"), 
           enrichment = factor(enrichment, c("enriched", "ns")))
  
# Significantly enriched GO terms ----------
  
  insert_msg("Significantly enriched GO terms")
  
  mibc_shared$go_significant <- mibc_shared$go_test %>% 
    filter(enrichment == "enriched") %>% 
    blast(pair_id) %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$term)
  
# Weighted similarity graph ---------
  
  insert_msg("Weighted similarity graph")

  ## data frame defining the vertices, edges, and edge attributes
  ## the edges are weighted with J similarity coefficients and take 
  ## the numbers of up- and downregulated genes as attributes
  
  mibc_shared$graph_data$edges <- mibc_shared$simil_test %>% 
    mutate(weight = j + 1e-6)
  
  mibc_shared$graph_data$gene_numbers <- mibc_shared$numbers %>% 
    select(pair_id, regulation, n) %>% 
    pivot_wider(id_cols = pair_id, 
                names_from = regulation, 
                values_from = n)
  
  mibc_shared$graph_data <- mibc_shared$graph_data %>% 
    reduce(left_join, by = "pair_id") %>% 
    mutate(upregulated = ifelse(is.na(upregulated), 0, upregulated), 
           downregulated = ifelse(is.na(downregulated), 0, downregulated), 
           edge_lab = paste("J =", signif(j, 2)), 
           edge_lab = paste(edge_lab, upregulated, sep = "\nup: n = "), 
           edge_lab = paste(edge_lab, downregulated, sep = "\ndown: n = "), 
           edge_lab = ifelse(j == 0, NA, edge_lab))

  ## the graph, setting the vertex attributes
  
  mibc_shared$simil_graph <- mibc_shared$graph_data %>% 
    graph_from_data_frame(directed = FALSE) %>% 
    set_vertex_attributes(mibc_shared$attr_df)
  
# Caching --------
  
  insert_msg("Caching")
  
  mibc_shared <- 
    mibc_shared[c("simil_test", "shared_features", "numbers", 
                   "go_test", "go_significant", "simil_graph")]
  
  save(mibc_shared, file = "./cache/mibc_shared.RData")
  
# END --------
  
  plan("sequential")
  
  insert_tail()