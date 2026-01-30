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
  
  ## numbers of common regulated genes 
  
  mibc_shared$total_numbers[c("bc", "consensus")] <- 
    list(mibc_bc, mibc_cons) %>% 
    map(~.x$numbers) %>% 
    map(filter, cohort == "common") %>% 
    map2(., c("n_clust_id", "n_consensusClass"), 
         ~mutate(.x, 
                 !!.y := n)) %>% 
    map(select, 
        regulation, 
        any_of(c("clust_id", "consensusClass", 
                 "n_clust_id", "n_consensusClass")))
  
  ## all and top 100 markers of bladder cancer clusters and consensus classes
  ## shared by at least three cohorts
  
  mibc_shared$markers <- 
    list(bc = mibc_bc, consensus = mibc_cons) %>% 
    map(~.x$cmm_markers) %>% 
    unlist(recursive = FALSE)
  
  mibc_shared$top_markers <- 
    list(bc = mibc_bc, consensus = mibc_cons) %>% 
    map(~.x$cmm_top_markers) %>% 
    unlist(recursive = FALSE)
  
  ## non-significant differences in consensus classes: 
  ## to be used for selection of distinct markers of the bladder cancer 
  ## clusters
  
  mibc_shared$anova_ns <- mibc_cons$anova %>% 
    map(filter, regulation == "ns") %>% 
    map(~.x$variable) %>% 
    shared_features(m = 4) %>% 
    as.character

  ## attributes for vertices of similarity graphs
  
  mibc_shared$attr_df <- globals$attr_df
  
  ## gene universe, to be used in GO enrichment analyses
  
  mibc_shared$universe <- mibc_bc$anova %>% 
    map(~.x$entrez_id) %>% 
    reduce(union)
  
  ## pairs of molecular subsets
  
  mibc_shared$pairs <- paste0("consensus.", globals$consensus_levels) %>% 
    map(function(x) paste0("bc.", globals$cluster_levels) %>% 
          map(~c(.x, x))) %>% 
    unlist(recursive = FALSE)
  
  ## input data: 
  ## vectors of simplified gene regulation information: 
  ## the upregulated and downregulated features are appended with suffixes
  
  mibc_shared$simil_data <- mibc_shared$genes %>% 
    map(map, ~map2(.x, names(.x), paste, sep = "|")) %>% 
    map(map, reduce, union) %>% 
    unlist(recursive = FALSE)
  
# Jaccard distances --------
  
  insert_msg("Jaccard distances")

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
                                   globals$consensus_levels))

# Numbers of genes in the intersections ---------
  
  insert_msg("Numbers of genes in the intersections")
  
  ## intersections
  
  mibc_shared$numbers <- mibc_shared$shared_features %>% 
    count(pair_id, regulation) %>% 
    left_join(mibc_shared$pair_df, by = "pair_id")
  
  ## appending with numbers of all genes regulated in particular 
  ## molecular subsets
  
  mibc_shared$numbers <- c(mibc_shared["numbers"], 
                           mibc_shared$total_numbers) %>% 
    reduce(left_join)

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
  
# labels for the overlaps between the subset schemes based on the GO enrichment -------
  
  insert_head("Labels for overlaps")
  
  mibc_shared$shared_labels <- 
    list(`bc.#1|consensus.Ba/Sq.downregulated` = c("biological process"), 
         `bc.#1|consensus.Ba/Sq.upregulated` = c("B & myeloid leukocytes", 
                                                 "phagocytosis, chemotaxis"), 
         `bc.#1|consensus.LumNS.downregulated` = c("epidermis, keratinocytes", 
                                                   "proteolysis"), 
         `bc.#1|consensus.LumNS.upregulated` = c("bones, muscles", 
                                                 "angiogenesis", 
                                                 "ECM, collagens", 
                                                 "MAPK, PI3K, IL1"), 
         `bc.#1|consensus.LumP.downregulated` = c("epidermis, keratinocytes"), 
         `bc.#1|consensus.LumU.downregulated` = c("epidermis, keratinocytes", 
                                                  "proteolysis"), 
         `bc.#1|consensus.NE-like.downregulated` = c("epithelium", 
                                                     "steroids"), 
         `bc.#1|consensus.Stroma-rich.downregulated` = c("epidermis, keratinocytes", 
                                                         "epithelium, adhesion"), 
         `bc.#1|consensus.Stroma-rich.upregulated` = c("EMT, organ development", 
                                                       "angiogenesis, migration", 
                                                       "ECM, collagens, integrins", 
                                                       "eicosanoids, retinoids", 
                                                       "TGFB, BMP, WNT, IL1", 
                                                       "PDGF, FGFR, MAPK, PI3K,", 
                                                       "small GTPases, TNF", 
                                                       "PLC, Ca2+ signaling"), 
         `bc.#2|consensus.Ba/Sq.downregulated` = c("EMT, organ development", 
                                                   "angiogenesis", 
                                                   "FAOX, lipids, steroids", 
                                                   "amines, retinoids, sugars", 
                                                   "TGFB, BMP, SMAD", 
                                                   "WNT, small GTPases"), 
         `bc.#2|consensus.Ba/Sq.upregulated` = c("epidermis, keratinocytes", 
                                                 "ECM, membranes", 
                                                 "DNA replication & damage", 
                                                 "cell cycle, cell death", 
                                                 "NK & T cells, leukocytes", 
                                                 "nucleotide metabolism", 
                                                 "glycolysis", 
                                                 "IFN/STAT, IL, TLR, NFkB", 
                                                 "EGFR/ERBB, MAPK"), 
         `bc.#2|consensus.NE-like.downregulated` = c("lipids, retinoids"),
         `bc.#2|consensus.NE-like.upregulated` = c("mitosis, cell cycle", 
                                                   "DNA replication & damage", 
                                                   "cell death"), 
         `bc.#2|consensus.Stroma-rich.downregulated` = c("vesicle transport"), 
         `bc.#2|consensus.Stroma-rich.upregulated` = c("organ development", 
                                                       "fibroblasts, migration", 
                                                       "ECM, collagens", 
                                                       "WNT, TGFB, BMP, SMAD", 
                                                       "MAPK, PI3K, small GTPases"), 
         `bc.#3|consensus.LumNS.downregulated` = c("epidermis, keratinocytes", 
                                                   "organ development", 
                                                   "cytokines, proteolysis", 
                                                   "ECM, collagens"), 
         `bc.#3|consensus.LumNS.upregulated` = c("lipids, steroids", 
                                                 "ion homeostasis"), 
         `bc.#3|consensus.LumP.downregulated` = c("epidermis, keratinocytes", 
                                                  "B, T & NK cells, leukocytes", 
                                                  "EMT, organs, angiogenesis", 
                                                  "chemokines, migration", 
                                                  "ECM, collagens, integrins", 
                                                  "IFN/STAT, IL, TNF, TLR", 
                                                  "NFkB, FcR, MAPK, PI3K", 
                                                  "WNT, BMP, PDGFR, VEGF", 
                                                  "GPCR, small GTPases"), 
         `bc.#3|consensus.LumP.upregulated` = c("epithelial development", 
                                                "fatty acid metabolism", 
                                                "FAOX, lipid metabolism", 
                                                "eicosanoids, steroids", 
                                                "retinoids, xenobiotics", 
                                                "amines, amino acids", 
                                                "ERBB, MAPK, small GTPases", 
                                                "WNT, NOTCH, LIF, BMP"), 
         `bc.#3|consensus.LumU.downregulated` = c("epidermis, keratinocytes", 
                                                  "ECM, MAPK"), 
         `bc.#3|consensus.LumU.upregulated` = c("biological process"), 
         `bc.#3|consensus.Stroma-rich.downregulated` = c("epidermis, keratinocytes", 
                                                         "proteolysis")) %>% 
    map_chr(paste, collapse = "\n")
  
# Identification of overlaps in cluster and consensus class markers ---------
  
  insert_msg("Identification of overlaps in markers")
  
  ## shared markers
  
  mibc_shared$shared_markers <- mibc_shared$top_markers %>% 
    find_overlaps(pairs = mibc_shared$pairs, as_list = TRUE) %>% 
    map(function(x) if(length(x) == 0) NULL else x) %>% 
    compact
  
  ## distinct markers for bladder cancer clusters #1, #2, #3, and, 
  ## respectively, stroma-rich, basal and luminal consensus class.
  ## additional selection step for the distinct markers: 
  ## lacking significance in ANOVA
  
  mibc_shared$distinct_markers <- 
    list(`#1` = c("bc.#1", "consensus.Stroma-rich"), 
         `#2` = c("bc.#2", "consensus.Ba/Sq"), 
         `#3` = c("bc.#3", "consensus.LumP")) %>% 
    map(~setdiff(mibc_shared$markers[[.x[[1]]]], 
                 mibc_shared$markers[[.x[[2]]]])) %>% 
    map(intersect, mibc_shared$anova_ns)

# Weighted similarity graph ---------
  
  insert_msg("Weighted similarity graph")

  ## data frame defining the vertices, edges, and edge attributes
  ## the edges are weighted with J similarity coefficients and take 
  ## the numbers of up- and downregulated genes, 
  ## labels of the gene overlaps with GO enrichment summaries, 
  ## and overlapping marker genes as attributes
  
  ## edge definitions and weights
  
  mibc_shared$graph_data$edges <- mibc_shared$simil_test %>% 
    mutate(weight = j + 1e-6)
  
  ## gene number attributes
  
  mibc_shared$graph_data$gene_numbers <- mibc_shared$numbers %>% 
    select(pair_id, regulation, n) %>% 
    pivot_wider(id_cols = pair_id, 
                names_from = regulation, 
                values_from = n)
  
  ## overlap labels
  
  mibc_shared$graph_data$labels <- mibc_shared$shared_labels %>% 
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
  
  mibc_shared$graph_data$markers <- mibc_shared$shared_markers %>% 
    map(sort) %>% 
    map(wrap_vector, len = 4) %>% 
    map_chr(paste, collapse = ", ") %>% 
    compress(names_to = "pair_id", 
             values_to = "marker_label")
  
  ## the whole graph definition
  
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
                  "go_test", "go_significant", "shared_labels", 
                  "shared_markers", "distinct_markers", 
                  "simil_graph")]
  
  save(mibc_shared, file = "./cache/mibc_shared.RData")
  
# END --------
  
  plan("sequential")
  
  insert_tail()