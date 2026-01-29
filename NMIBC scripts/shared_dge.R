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
# Significant enrichment considered for at least 3 genes assigned to the term, 
# pFDR < 0.05, and OR < 0.05.
#
# 4) Construction of a weighted similarity graph with Jaccard similarities.
#
# 5) Labels for the ovelaps with key biological processes derived from 
# the GO enrichment analyses

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
  
  ## numbers of common regulated genes 
  
  nmibc_shared$total_numbers[c("bc", "uromol")] <- 
    list(nmibc_bc, nmibc_uro) %>% 
    map(~.x$numbers) %>% 
    map(filter, cohort == "common") %>% 
    map2(., c("n_clust_id", "n_NMIBC_class"), 
         ~mutate(.x, 
                 !!.y := n)) %>% 
    map(select, 
        regulation, 
        any_of(c("clust_id", "NMIBC_class", 
                 "n_clust_id", "n_NMIBC_class")))
  
  ## common top markers of bladder cancer clusters and UROMOL classes
  
  nmibc_shared$markers <- 
    list(bc = nmibc_bc, uromol = nmibc_uro) %>% 
    map(~.x$cmm_top_markers) %>% 
    unlist(recursive = FALSE)
  
  ## attributes for vertices of similarity graphs
  
  nmibc_shared$attr_df <- globals$attr_df
  
  ## gene universe, to be used in GO enrichment analyses
  
  nmibc_shared$universe <- nmibc_bc$anova %>% 
    map(~.x$entrez_id) %>% 
    reduce(union)
  
  ## pairs of molecular subsets
  
  nmibc_shared$pairs <- paste0("uromol.", globals$uromol_levels) %>% 
    map(function(x) paste0("bc.", globals$cluster_levels) %>% 
          map(~c(.x, x))) %>% 
    unlist(recursive = FALSE)
  
  ## input data for similarity computation: 
  ## vectors of simplified gene regulation information: 
  ## the upregulated and downregulated features are appended with suffixes
  
  nmibc_shared$simil_data <- nmibc_shared$genes %>% 
    map(map, ~map2(.x, names(.x), paste, sep = "|")) %>% 
    map(map, reduce, union) %>% 
    unlist(recursive = FALSE)
  
# Jaccard distances --------
  
  insert_msg("Jaccard distances")

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
  
  ## appending with numbers of all genes regulated in particular 
  ## molecular subsets
  
  nmibc_shared$numbers <- c(nmibc_shared["numbers"], 
                            nmibc_shared$total_numbers) %>% 
    reduce(left_join)
  
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
  
# labels for the overlaps between the subset schemes based on the GO enrichment -------
  
  insert_head("Labels for overlaps")
  
  nmibc_shared$shared_labels <- 
    list(`bc.#1|uromol.2b.upregulated` = c("antigens, IL, TNF", 
                                           "B, T, NK cells, leukocytes", 
                                           "membrane organization", 
                                           "migration, phagocytosis",
                                           "angiogenesis", 
                                           "cell death, ECM, integrins", 
                                           "IFN/STAT, FcR, PI3K", 
                                           "VEGF, small GTPases", 
                                           "GPCR, PLC", 
                                           "WNT, MAPK, FGFR, NFkB", 
                                           "eicosanoids, Ca2+"), 
         `bc.#2|uromol.2b.upregulated` = c("antigens, IL", 
                                           "B, T, NK cells & leukocytes", 
                                           "migration, phagocytosis", 
                                           "cell death, angiogenesis", 
                                           "IFN/STAT, FcR, NFkB", 
                                           "GPCR, PLC, MAPK", 
                                           "small GTPases"), 
         `bc.#2|uromol.2b.downregulated` = c("miRNA"), 
         `bc.#3|uromol.1.downregulated` = c("bone & organ devlopment", 
                                            "epithelial differentiation", 
                                            "adhesion, migration, EMT", 
                                            "angiogenesis, cell death", 
                                            "WNT, BMP, chemokines", 
                                            "MAPK, small GTPases", 
                                            "TNF, IL"), 
         `bc.#3|uromol.2a.downregulated` = c("B, T, NK cells", 
                                             "cell death", 
                                             "ECM, collagens", 
                                             "IFN, IL, TLR, TNF", 
                                             "small GTPases, MAPK", 
                                             "Ca2+ signaling"), 
         `bc.#3|uromol.3.downregulated` = c("antigens, lectins", 
                                            "B, T, NK cells", 
                                            "chemotaxis, membrane", 
                                            "angiogenesis, cell death",
                                            "ECM, collagens", 
                                            "eicosanoids, vit. D", 
                                            "FcR, IFN/STAT, IL, TLR", 
                                            "WNT, NFkB, MAPK, PI3K", 
                                            "GPCR, PLC, small GTPases", 
                                            "Ca2+ & synaptic signaling")) %>% 
    map_chr(paste, collapse = "\n")
  
# Identification of overlaps in cluster and UROMOL class markers ---------
  
  insert_msg("Identification of overlaps in markers")
  
  nmibc_shared$shared_markers <- nmibc_shared$markers %>% 
    find_overlaps(pairs = nmibc_shared$pairs, 
                  as_list = TRUE) %>% 
    map(function(x) if(length(x) == 0) NULL else x) %>% 
    compact
  
# Weighted similarity graph ---------
  
  insert_msg("Weighted similarity graph")

  ## data frame defining the vertices, edges, and edge attributes
  ## the edges are weighted with J similarity coefficients and take 
  ## the numbers of up- and downregulated genes, 
  ## labels of the gene overlaps with GO enrichment summaries, 
  ## and text with overlapping markers as attributes
  
  ## edge definitions and weights
  
  nmibc_shared$graph_data$edges <- nmibc_shared$simil_test %>% 
    mutate(weight = j + 1e-6)
  
  ## gene number attributes of the edges
  
  nmibc_shared$graph_data$gene_numbers <- nmibc_shared$numbers %>% 
    select(pair_id, regulation, n) %>% 
    pivot_wider(id_cols = pair_id, 
                names_from = regulation, 
                values_from = n)
  
  ## overlap labels
  
  nmibc_shared$graph_data$labels <- nmibc_shared$shared_labels %>% 
    compress(names_to = "condition", 
             values_to = "label") %>% 
    mutate(regulation = stri_extract(condition, 
                                     regex = "upregulated|downregulated"), 
           regulation = paste0("label_", regulation), 
           pair_id = stri_replace(condition, 
                                  regex = "\\.(upregulated|downregulated)", 
                                  replacement = "")) %>% 
    select(-condition) %>% 
    pivot_wider(id_cols = pair_id, 
                names_from = regulation, 
                values_from = label)
  
  ## overlapping markers
  
  nmibc_shared$graph_data$markers <- nmibc_shared$shared_markers %>% 
    map(sort) %>% 
    map(wrap_vector) %>% 
    map_chr(paste, collapse = ", ") %>% 
    compress(names_to = "pair_id", 
             values_to = "marker_label")
  
  ## the whole definition
  
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
                   "go_test", "go_significant", 
                   "shared_labels", "shared_markers", 
                   "simil_graph")]
  
  save(nmibc_shared, file = "./cache/nmibc_shared.RData")
  
# END --------
  
  plan("sequential")
  
  insert_tail()