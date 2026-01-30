# Tables with summaries of results of differential gene expression analyses, 
# marker search, and GO enrichment analyses for the overlapping genes 

  insert_head()
  
# container --------
  
  sum_tabs <- list()
  tab_globals <- list()
  
# table globals --------
  
  insert_msg("Table globals")
  
  tab_globals$systems <- 
    c("bladder cancer clusters, NMIBC", 
      "UROMOL classes, NMIBC" , 
      "bladder cancer clusters, MIBC",  
      "consensus classes, MIBC")
  
  ## data frames with pair IDs and molecular subsets
  ## used later in formatting of GO enrichment and gene counting results
  
  tab_globals$pair_df <- 
    list(nmibc = nmibc_shared, 
         mibc = mibc_shared) %>% 
    map(~.x$simil_test) %>% 
    map(select, 
        pair_id, 
        clust_id, 
        any_of(c("NMIBC_class", "consensusClass"))) %>% 
    map(filter, !duplicated(pair_id))
  
# differential gene expression --------
  
  insert_msg("Differential gene expression")
  
  sum_tabs$dge_systems <- 
    list(nmibc_bc, nmibc_uro, mibc_bc, mibc_cons) %>% 
    set_names(tab_globals$systems) %>% 
    map(function(tst) tst$test %>% 
          map(filter, 
              variable %in% unique(unlist(tst$cmm_significant))) %>% 
          map(mutate, entrez_id = as.character(entrez_id)) %>% 
          map2(map(tst$anova, 
                   transmute, 
                   variable = variable, 
                   p_anova = p_adjusted, 
                   etasq = etasq), 
               left_join, 
               by = "variable") %>% 
          compress(names_to = "cohort"))
  
  sum_tabs$dge_systems <- 
    map2(sum_tabs$dge_systems, 
         c("clust_id", "NMIBC_class", "clust_id", "consensusClass"), 
         ~mutate(.x, 
                 subset = .data[[.y]])) %>% 
    map(~.x[!names(.x) %in% c("clust_id", 
                              "NMIBC_class", 
                              "clust_id", 
                              "consensusClass")]) %>% 
    compress(names_to = "system")
  
  sum_tabs$dge_systems <- sum_tabs$dge_systems %>% 
    transmute(`Molecular classification, cancer entity` = 
                factor(system, tab_globals$systems), 
              Cohort = globals$cohort_labs[cohort], 
              `Gene symbol` = gene_symbol, 
              `Entrez ID` = entrez_id, 
              `Molecular subset` = subset, 
              `Regulation status vs cohort mean` = regulation, 
              `ANOVA, pFDR` = signif(p_anova, 2), 
              `ANOVA, effect size, eta-square` = signif(etasq, 2), 
              `log2 FC regulation vs cohort mean, 95% CI` = 
                paste0(signif(deviation_center, 2), 
                      " [", signif(lower_ci, 2), " to ", 
                      signif(upper_ci, 2), "]"), 
              `post-hoc T test, pFDR` = signif(p_adjusted, 2)) %>% 
    arrange(`Molecular classification, cancer entity`, 
            Cohort, 
            `Gene symbol`)
  
  sum_tabs$dge_systems <-  sum_tabs$dge_systems %>% 
    as_mdtable(label = "dge_systems", 
               ref_name = "dge_systems", 
               caption = paste("Differential gene expression in bladder cancer", 
                               "clusters of NMIBC, UROMOL molecular classes", 
                               "of NMIBC, and consensus classes of MIBC.", 
                               "The differential gene expression in the molecular", 
                               "subsets as compared with the cohort average was", 
                               "investigated by one-way ANOVA with eta-square", 
                               "effect size statistic and one-sample post-hoc", 
                               "T test.", 
                               "P values were corrected for multiple testing", 
                               "with the false discovery rate (FDR) method.", 
                               "Significant regulation was considered for", 
                               "ANOVA pFDR < 0.05, eta-square >= 0.06, and", 
                               "T test pFDR < 0.05.", 
                               "Log2 fold-change (FC) estimates with 95% confidence", 
                               "intervals (95%CI), and FDR-adjusted p values in", 
                               "one-way ANOVA and post-hoc tests are listed for",
                               "genes found to be significantly", 
                               "regulated in at least three cohorts.", 
                               "The table is available as a supplementary", 
                               "Excel file."))
  
# differentially regulated genes shared by the subtype schemes -------
  
  insert_msg("Regulated genes shared by the subtypes")
  
  sum_tabs$shared_genes_systems <- 
    list(NMIBC = nmibc_shared, MIBC = mibc_shared) %>% 
    map(~.x$shared_features) %>% 
    map2(., c("NMIBC_class", "consensusClass"),
         ~transmute(.x, 
                    `Bladder cancer cluster` = clust_id, 
                    `UROMOL/consensus class` = .data[[.y]], 
                    `Gene symbol` = gene_symbol, 
                    `Entrez ID` = entrez_id, 
                    `Regulation status vs cohort mean` = regulation)) %>% 
    compress(names_to = "Cancer entity") %>% 
    relocate(`Cancer entity`)
  
  sum_tabs$shared_genes_systems <- sum_tabs$shared_genes_systems %>% 
    as_mdtable(label = "shared_genes_systems", 
               ref_name = "shared_genes_systems", 
               caption = paste("Genes differentially regulated in at least", 
                               "three cohorts shared by", 
                               "bladder cancer clusters and UROMOL classes", 
                               "in NMIBC, and bladder cancer clusters and", 
                               "consensus classes in NMIBC.",
                               "Molecular subset names, and regulation status", 
                               "in the molecular subsets as compared with the", 
                               "cohort mean are listed.", 
                               "The table is available as a supplementary", 
                               "Excel file."))
  
# Jaccard similarity coefficients --------
  
  insert_msg("Jaccard similarity coefficients")
  
  sum_tabs$jaccard_systems <- 
    list(NMIBC = nmibc_shared,
         MIBC = mibc_shared) %>% 
    map2(., tab_globals$pair_df, 
         ~list(.y, 
               .x$numbers %>% 
                 pivot_wider(id_cols = pair_id, 
                             names_from = regulation, 
                             values_from = n), 
               .x$simil_test[, c("pair_id", "j")])) %>% 
    map(reduce, left_join, by = "pair_id")
  
  sum_tabs$jaccard_systems <- 
    map2(sum_tabs$jaccard_systems, 
         c("NMIBC_class", "consensusClass"), 
         ~transmute(.x, 
                    `Bladder cancer cluster` = clust_id, 
                    `UROMOL/consensus class` = .data[[.y]], 
                    `Shared upregulated genes, N` = ifelse(is.na(upregulated), 
                                                           0, upregulated), 
                    `Shared downregulated genes, N` = ifelse(is.na(downregulated), 
                                                             0, downregulated), 
                    `Similarity coefficient, J` = signif(j, 2))) %>% 
    compress(names_to = "Cancer entity")
  
  sum_tabs$jaccard_systems <- sum_tabs$jaccard_systems %>% 
    as_mdtable(label = "jaccard_systems", 
               ref_name = "jaccard_systems", 
               caption = paste("Similarity of bladder cancer clusters and", 
                               "UROMOL classes of NMIBC, bladder cancer clusters", 
                               "and consensus classes of NMIBC in respect to", 
                               "shared differentially regulated genes.", 
                               "Numbers of concordantly up- and downregulated genes", 
                               "in the molecular subsets as compared with the", 
                               "cohort average, and Jaccard similarity coefficients", 
                               "are listed."))
  
# GO enrichment analyses for the overlapping genes --------
  
  insert_msg("GO enrichment analyses for the overlapping regulated genes")
  
  sum_tabs$go_systems <- 
    list(NMIBC = nmibc_shared, 
         MIBC = mibc_shared) %>% 
    map(~.x$go_test) %>% 
    map(filter, enrichment == "enriched")
  
  sum_tabs$go_systems <- 
    map2(sum_tabs$go_systems, 
         c("NMIBC_class", "consensusClass"), 
         ~transmute(.x, 
                    `Bladder cancer cluster` = clust_id, 
                    `UROMOL/consensus class` = .data[[.y]], 
                    `Regulation status vs cohort mean` = regulation, 
                    `GO term name` = term, 
                    `GO term ID` = go_id, 
                    `Differentially regulated genes in term, N` = n_go_dge, 
                    `Odds ratio` = signif(or, 2), 
                    `pFDR` = signif(p_adjusted, 2))) %>% 
    compress(names_to = "Cancer entity") %>% 
    relocate(`Cancer entity`)
  
  sum_tabs$go_systems <- sum_tabs$go_systems %>% 
    as_mdtable(label = "go_systems", 
               ref_name = "go_systems", 
               caption = paste("Biological process gene ontology (GO) enrichment", 
                               "analysis for genes concordantly regulated in", 
                               "bladder cancer clusters and UROMOL classes of NMIBC,", 
                               "and bladder cancer clusters and consensus", 
                               "classes of MIBC.", 
                               "Enrichment magnitude in the shared gene sets as", 
                               "compared with the entire genome was measured by", 
                               "the odds ratio statistic.", 
                               "P values were corrected for multiple testing with", 
                               "the false discovery rate (FDR) method.", 
                               "Significant enrichment was considered for odds", 
                               "ratio >= 1.44, at least three regulated genes", 
                               "per term, and pFDR < 0.05.", 
                               "Results for significantly enriched GO terms are", 
                               "listed.", 
                               "The table is available as a supplementary", 
                               "Excel file."))
  
# marker search -----------
  
  insert_msg("Marker search")
  
  sum_tabs$marker_systems <- 
    list(nmibc_bc, nmibc_uro, mibc_bc, mibc_cons) %>% 
    set_names(tab_globals$systems) %>% 
    map(function(tst) map(tst$roc, 
                          filter, 
                          variable %in% unique(unlist(tst$cmm_markers)))) %>% 
    map(map, mutate, entrez_id = as.character(entrez_id)) %>% 
    map(compress, names_to = "cohort")
  
  sum_tabs$marker_systems <- 
    map2(sum_tabs$marker_systems, 
         c("clust_id", "NMIBC_class", "clust_id", "consensusClass"), 
         ~mutate(.x, 
                 subset = .data[[.y]])) %>% 
    map(~.x[!names(.x) %in% c("clust_id", 
                              "NMIBC_class", 
                              "clust_id", 
                              "consensusClass")]) %>% 
    compress(names_to = "system")
  
  sum_tabs$marker_systems <- sum_tabs$marker_systems %>% 
    transmute(`Molecular classification, cancer entity` = 
                factor(system, tab_globals$systems), 
              `Molecular subset` = subset, 
              Cohort = globals$cohort_labs[cohort],
              `Gene symbol` = gene_symbol, 
              `Entrez ID` = entrez_id, 
              `Marker, AUC >= 0.714` = marker, 
              `log2 expression cutoff` = signif(cutoff, 3), 
              AUC = signif(auc, 2), 
              Sensitivity = signif(Se, 2), 
              Specificity = signif(Sp, 2)) %>% 
    arrange(`Molecular classification, cancer entity`, 
            `Molecular subset`, 
            `Gene symbol`, 
            Cohort)
  
  sum_tabs$marker_systems <- sum_tabs$marker_systems %>% 
    as_mdtable(label = "marker_systems", 
               ref_name = "marker_systems", 
               caption = paste("Marker genes of bladder cancer clusters of NMIBC,", 
                               "UROMOL classes of NMIBC, bladder cancer clusters", 
                               "of MIBC, and consensus classes of MIBC.", 
                               "The marker genes were identified by false", 
                               "discovery rate (FDR) corrected one-way ANOVA", 
                               "with eta-square effect size statistic, and", 
                               "receiver-operating characteristic.", 
                               "Markers of the molecular subsets were identified", 
                               "by significant differences in expression", 
                               "with ANOVA pFDR < 0.05 and eta-square >= 0.06,", 
                               "and area under the ROC curve (AUC) >= 0.714.", 
                               "The optimal cutoffs of log2 gene expression levels", 
                               "for identification of the molecular subset", 
                               "found with the Youden criterion, sensitivity,", 
                               "specificity, and AUC are listed for", 
                               "the marker genes shared by at least three cohorts.", 
                               "The table is available as a supplementary Excel", 
                               "file."))

# Saving the tables on the disc --------
  
  insert_msg("Saving the tables on the disc")
  
  sum_tabs %>% 
    save_excel(path = "./report/tables.xlsx")
  
# END --------
  
  rm(tab_globals)
  
  insert_tail()