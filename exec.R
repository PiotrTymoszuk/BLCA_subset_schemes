# Launches the entire pipeline

  library(soucer)
  
  print(source_all(c("import.R", 
                     "NMIBC.R"), 
                   message = TRUE, crash = TRUE))