# Launches the entire pipeline

  library(soucer)
  
  print(source_all(c("import.R", 
                     "NMIBC.R", 
                     "MIBC.R", 
                     "summary.R"), 
                   message = TRUE, crash = TRUE))