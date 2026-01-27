# Similarity and dissimilarity of the bladder cancer clusters and UROMOL 
# classes in NMIBC: 
#
# 1) Jaccard similarities between the clusters and UROMOL classes in respect 
# to genes found differentially regulated in at least three cohorts.
#
# 2) Identification of shared upregulated and downregulated genes 
# between the bladder cancer clusters and UROMOL classes
#
# 3) GO enrichment analyses for the overlapping up- and downregulated genes. 
# Significant enrichment considere for at least 3 genes assigned to the term, 
# pFDR < 0.05, and OR < 0.05.
#
# 4) Construction of a weighted similarity graph with Jaccard similarities.

  insert_head()
  
# container -------
  
  nmibc_shared <- list()
  
# parallel backend --------
  
  insert_msg("Parallel backend")
  
  plan("multisession")
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  ## common differentially regulated genes in bladder cancer clusters 
  ## and common regulated genes in UROMOL classes
  
  nmibc_shared$genes[c("bc", "uromol")] <- 
    list(nmibc_bc, nmibc_uro) %>% 
    map(~.x$cmm_significant)
  
  ## attributes for vertices of similarity graphs
  
  nmibc_shared$attr_df <- globals$attr_df
  
  ## gene universe, to be used in GO enrichment analyses
  
  nmibc_shared$universe <- nmibc_bc$anova %>% 
    map(~.x$entrez_id) %>% 
    reduce(union)
  
# Jaccard distances --------
  
  insert_msg("Jaccard distances")
  
  ## input data: vectors of simplified gene regulation information: 
  ## the upregulate and downregulated features are appended with suffixes
  
  nmibc_shared$simil_data <- nmibc_shared$genes %>% 
    map(map, ~map2(.x, names(.x), paste, sep = "|")) %>% 
    map(map, reduce, union) %>% 
    unlist(recursive = FALSE)
  
  ## Jaccard distances
  
  nmibc_shared$simil_test <- nmibc_shared$simil_data %>% 
    set_similarity(method = "jaccard")
  
  nmibc_shared$simil_test <- 
    nmibc_shared$simil_test[c("bc.#1", "bc.#2", "bc.#3"), 
                            c("uromol.1", "uromol.2a", "uromol.2b", "uromol.3")] %>% 
    mtx2long(row_var = "clust_id", 
             col_var = "NMIBC_class", 
             value_var = "j") %>% 
    mutate(clust_id = stri_split_fixed(clust_id, 
                                       pattern = ".", 
                                       simplify = TRUE)[, 2], 
           clust_id = factor(clust_id, globals$cluster_levels), 
           NMIBC_class = stri_split_fixed(NMIBC_class, 
                                          pattern = ".", 
                                          simplify = TRUE)[, 2], 
           NMIBC_class = factor(NMIBC_class, globals$uromol_levels), 
           pair_id = paste(paste0("bc.", clust_id), 
                           paste0("uromol.", NMIBC_class), 
                           sep = "|"))
  
  ## a data frame with pair IDs and molecular subsets
  ## used later in formatting of GO enrichment and gene counting results

  nmibc_shared$pair_df <- 
    nmibc_shared$simil_test[, c("pair_id", "clust_id", "NMIBC_class")] %>% 
    filter(!duplicated(pair_id))

# Identification of shared differentially regulated genes ---------
  
  insert_msg("Identification of shared regulated genes")
  
  ## pairs of molecular subsets
  
  nmibc_shared$pairs <- paste0("uromol.", globals$uromol_levels) %>% 
    map(function(x) paste0("bc.", globals$cluster_levels) %>% 
          map(~c(.x, x))) %>% 
    unlist(recursive = FALSE)
  
  ## common regulated features, appended with Entrez IDs

  nmibc_shared$shared_features <- 
    find_overlaps(nmibc_shared$simil_data, 
                  nmibc_shared$pairs) %>% 
    mutate(gene_symbol = variable, 
           entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              keytype = "SYMBOL", 
                              column = "ENTREZID"), 
           clust_id = stri_split_fixed(element1, 
                                       pattern = ".", 
                                       simplify = TRUE)[, 2], 
           clust_id = factor(clust_id, globals$cluster_levels), 
           NMIBC_class = stri_split_fixed(element2, 
                                          pattern = ".", 
                                          simplify = TRUE)[, 2], 
           NMIBC_class = factor(NMIBC_class, globals$uromol_levels))

# Numbers of genes in the intersections ---------
  
  insert_msg("Numbers of genes in the intersections")
  
  nmibc_shared$numbers <- nmibc_shared$shared_features %>% 
    count(pair_id, regulation) %>% 
    left_join(nmibc_shared$pair_df, by = "pair_id")
  
# GO enrichment analyses for the overlaps with at least 10 genes ----------
  
  insert_msg("GO enrichment analyses for the overlaps")
  
  nmibc_shared$go_test <- nmibc_shared$shared_features %>% 
    blast(pair_id, regulation) %>% 
    map(~.x$entrez_id) %>% 
    map(function(x) if(length(x) < 10) NULL else x) %>% 
    compact %>% 
    future_map(GOana, 
               universe = nmibc_shared$universe, 
               ontology = "BP", 
               .options = furrr_options(seed = 1232))
  
  ## formatting the GO enrichment results
  
  nmibc_shared$go_test <- nmibc_shared$go_test %>% 
    compress(names_to = "condition") %>% 
    mutate(regulation = stri_extract(condition, 
                                     regex = "upregulated|downregulated"), 
           regulation = factor(regulation, globals$regulation_levels), 
           pair_id = stri_replace(condition, 
                                  regex = "\\.(upregulated|downregulated)", 
                                  replacement = "")) %>% 
    left_join(nmibc_shared$pair_df, 
              by = "pair_id")
  
  ## identification of significant GO terms
  
  nmibc_shared$go_test <- nmibc_shared$go_test %>% 
    filter(n_go_dge >= 3) %>% 
    mutate(enrichment = ifelse(p_adjusted >= 0.05 | or < 1.44, 
                               "ns", "enriched"), 
           enrichment = factor(enrichment, c("enriched", "ns")))
  
# Significantly enriched GO terms ----------
  
  insert_msg("Significantly enriched GO terms")
  
  nmibc_shared$go_significant <- nmibc_shared$go_test %>% 
    filter(enrichment == "enriched") %>% 
    blast(pair_id) %>% 
    map(blast, regulation) %>% 
    map(map, ~.x$term)
  
# Weighted similarity graph ---------
  
  insert_msg("Weighted similarity graph")

  ## data frame defining the vertices, edges, and edge attributes
  ## the edges are weighted with J similarity coefficients and take 
  ## the numbers of up- and downregulated genes as attributes
  
  nmibc_shared$graph_data$edges <- nmibc_shared$simil_test %>% 
    mutate(weight = j + 1e-6)
  
  nmibc_shared$graph_data$gene_numbers <- nmibc_shared$numbers %>% 
    select(pair_id, regulation, n) %>% 
    pivot_wider(id_cols = pair_id, 
                names_from = regulation, 
                values_from = n)
  
  nmibc_shared$graph_data <- nmibc_shared$graph_data %>% 
    reduce(left_join, by = "pair_id") %>% 
    mutate(upregulated = ifelse(is.na(upregulated), 0, upregulated), 
           downregulated = ifelse(is.na(downregulated), 0, downregulated), 
           edge_lab = paste("J =", signif(j, 2)), 
           edge_lab = paste(edge_lab, upregulated, sep = "\nup: n = "), 
           edge_lab = paste(edge_lab, downregulated, sep = "\ndown: n = "), 
           edge_lab = ifelse(j == 0, NA, edge_lab))

  ## the graph, setting the vertex attributes
  
  nmibc_shared$simil_graph <- nmibc_shared$graph_data %>% 
    graph_from_data_frame(directed = FALSE) %>% 
    set_vertex_attributes(nmibc_shared$attr_df)
  
# Caching --------
  
  insert_msg("Caching")
  
  nmibc_shared <- 
    nmibc_shared[c("simil_test", "shared_features", "numbers", 
                   "go_test", "go_significant", "simil_graph")]
  
  save(nmibc_shared, file = "./cache/nmibc_shared.RData")
  
# END --------
  
  plan("sequential")
  
  insert_tail()