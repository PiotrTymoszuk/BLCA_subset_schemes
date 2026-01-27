# Analyses for bladder cancer clusters and consensus classes in NMIBC cohorts 
# (TCGA BLCA, IMvigor, GSE83586, GSE87304, GSE128192, GSE203149). 
#
# 1) Differential gene expression in bladder cancer clusters as compared with 
# the cohort average. One-way ANOVA and one-sample post-hoc T test. 
# Significant effects: pFDR ANOVA < 0.05, eta-square >= 0.06, pFDR T test < 0.05.
#
# 2) Differential gene expression in consensus classes as compared with 
# the cohort average. One-way ANOVA and one-sample post-hoc T test. 
# Significant effects: pFDR ANOVA < 0.05, eta-square >= 0.06, pFDR T test < 0.05.
#
# 3) Identification of differentially regulated gene shared by at least three 
# cohorts and between the molecular subset schemes. Quantification of the overlaps 
# in gene expression between the schemes with Jaccard's similarity coefficients. 
# GO enrichment analyses for the overlapping up- and downregulated genes.

# tools ---------
  
  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)
  
  library(microViz)
  library(fastTest)

  library(igraph)
  library(graphExtra)

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  library(furrr)
  
  library(soucer)

  reduce <- purrr::reduce
  select <- dplyr::select

  insert_head()
  
  c("./tools/globals.R", 
    "./tools/functions.R") %>% 
    source_all(message = TRUE, crash = TRUE)
  
# analysis globals --------
  
  insert_msg("Analysis globals")
  
  c("./MIBC scripts/globals.R") %>% 
    source_all(message = TRUE, crash = TRUE)

# analysis scripts --------
  
  insert_msg("Analysis scripts")
  
  ## cached results of differential expression analyses, and 
  ## similarity and dissimilarity of the bladder cancer clusters 
  ## and UROMOL classes
  
  list(cache_path = c("./cache/mibc_bc.RData", 
                      "./cache/mibc_cons.RData", 
                      "./cache/mibc_shared.RData"), 
       script_path = c("./MIBC scripts/bc_dge.R", 
                       "./MIBC scripts/consensus_dge.R", 
                       "./MIBC scripts/shared_dge.R"), 
       message = paste("Cached", 
                       c("DGE analysis results for bladder cancer clusters", 
                         "DGE analysis results for consensus classes", 
                         "overlaps of differentially regulated genes"))) %>% 
    pwalk(access_cache)
  
# END --------
  
  insert_tail()