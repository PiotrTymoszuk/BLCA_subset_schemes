# Import of cleared data for selected large cohorts of 
# NMIBC (non-muscle invasive bladder cancer) and 
# MIBC (muscle invasive bladder cancer): 
#
# 1) Assignment of NMIBC and MIBC samples in selected cohorts to bladder 
# cancer clusters (all samples), UROMOL classes (NMIBC only), and consensus 
# classes of MIBC (only MIBC)
#
# 2) The corresponding gene expression data.
#
# This information is imported from RData files.

# tools -------

  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)

  library(soucer)

  insert_head()
  
  c("./tools/globals.R", 
    "./tools/functions.R") %>% 
    source_all(message = TRUE, crash = TRUE)
  
# NMIBC and MIBC data from RData files ---------
  
  insert_msg("Reading the assignment and gene expression data")
  
  load("./data/ex_nmibc.RData")
  load("./data/ex_mibc.RData")
  
# END --------
  
  insert_tail()