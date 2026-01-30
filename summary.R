# Summary plots and tables to be exported to the main bladder cancer 
# cluster project

# tools ---------

  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)
  
  library(microViz)
  library(fastTest)
  library(exda)
  
  library(igraph)
  library(graphExtra)
  library(ggnetwork)

  library(ggtext)
  library(figur)
  library(cowplot)
  
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
  
# analysis scripts --------
  
  insert_msg("Analysis scripts")
  
  ## plots 
  
  c("./summary scripts/dge_plots.R", 
    "./summary scripts/graph_plots.R", 
    "./summary scripts/venn_plots.R", 
    "./summary scripts/text_plots.R") %>% 
    source_all(message = TRUE, crash = TRUE)
  
  ## figures and tables
  
  c("./summary scripts/figures.R", 
    "./summary scripts/tables.R") %>% 
    source_all(message = TRUE, crash = TRUE)
  
# exports --------
  
  insert_msg("Exports")
  
  save(sum_figs, file = "./exports/sum_figs.RData")
  save(sum_tabs, file = "./exports/sum_tabs.RData")
  
# END -------
  
  insert_tail()