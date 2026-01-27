# Globals for the project: graphics, labels, levels, colors, and shapes

# container ------

  globals <- list()

# graphic globals ------

  ## graphic theme for most of the plots

  globals$common_text <- element_text(size = 8,
                                      face = "plain",
                                      color = "black")
  
  globals$common_margin <- ggplot2::margin(t = 5,
                                           l = 4,
                                           r = 2,
                                           unit = "mm")
  
  globals$common_theme <- theme_classic() +
    theme(axis.text = globals$common_text,
          axis.title = globals$common_text,
          plot.title = element_text(size = 8,
                                    face = "bold"),
          plot.subtitle = globals$common_text,
          plot.tag = element_text(size = 8,
                                  face = "plain",
                                  color = "black",
                                  hjust = 0,
                                  vjust = 1),
          plot.tag.position = "bottom",
          legend.text = globals$common_text,
          legend.title = globals$common_text,
          strip.text = globals$common_text,
          strip.background = element_rect(fill = "gray95",
                                          color = "gray80"),
          plot.margin = globals$common_margin,
          panel.grid.major = element_line(color = "gray90"))

  ## graphic theme for network plots
  
  globals$net_theme <- theme_void() + 
    theme(plot.title = element_text(size = 8,
                                    face = "bold"), 
          plot.subtitle = globals$common_text, 
          legend.text = globals$common_text,
          legend.title = globals$common_text, 
          plot.margin = globals$common_margin)
  
# cohort labels ------

  globals$cohort_labs <- c(tcga = "TCGA BLCA",
                           imvigor = "IMvigor",
                           gse13507 = "GSE13507",
                           gse32548 = "GSE32548",
                           gse48075 = "GSE48075",
                           gse48276 = "GSE48276",
                           gse83586 = "GSE83586",
                           gse86411 = "GSE86411",
                           gse87304 = "GSE87304",
                           gse120736 = "GSE120736",
                           gse124305 = "GSE124305",
                           gse128192 = "GSE128192",
                           gse128701 = "GSE128701",
                           gse128959 = "GSE128959",
                           gse169455 = "GSE169455",
                           gse198269 = "GSE198269",
                           gse203149 = "GSE203149",
                           emtab4321 = "E-MTAB-4321")
  
# Large color palettes --------
  
  globals$tableau20_colors <- c(
    "Dark_Blue" = "#1F77B4",
    "Light_Blue" = "#AEC7E8",
    "Dark_Orange" = "#FF7F0E",
    "Light_Orange" = "#FFBB78",
    "Dark_Green" = "#2CA02C",
    "Light_Green" = "#98DF8A",
    "Dark_Red" = "#D62728",
    "Light_Red" = "#FF9896",
    "Dark_Purple" = "#9467BD",
    "Light_Purple" = "#C5B0D5",
    "Dark_Brown" = "#8C564B",
    "Light_Brown" = "#C49C94",
    "Dark_Pink" = "#E377C2",
    "Light_Pink" = "#F7B6D2",
    "Dark_Gray" = "#7F7F7F",
    "Light_Gray" = "#C7C7C7",
    "Dark_Yellow-Green" = "#BCBD22",
    "Light_Yellow-Green" = "#DBDB8D",
    "Dark_Teal" = "#17BECF",
    "Light_Teal" = "#9EDAE5"
  )
  
  globals$tableau10_colors <- c(
    "Blue" = "#4E79A7",
    "Orange" = "#F28E2B",
    "Red" = "#E15759",
    "Teal" = "#76B7B2",
    "Green" = "#59A14F",
    "Yellow" = "#EDC949",
    "Purple" = "#AF7AA1",
    "Pink" = "#FF9DA7",
    "Brown" = "#9C755F",
    "Gray" = "#BAB0AB"
  )
  
  globals$wes_grandbudapest1_colors <- c(
    "Pink" = "#F1BB7B",
    "Red" = "#FD6467",
    "Brown" = "#5B1A18",
    "Orange" = "#D67236"
  )
  
  globals$wes_grandbudapest2_colors <- c(
    "Light_Pink" = "#E6A0C4",
    "Light_Blue" = "#C6CDF7",
    "Tan" = "#D8A499",
    "Dark_Blue" = "#7294D4"
  )
  
  globals$wes_zissou1_colors <- c(
    "Light_Blue" = "#3B9AB2",
    "Medium_Blue" = "#78B7C5",
    "Yellow" = "#EBCC2A",
    "Dark_Yellow" = "#E1AF00",
    "Orange_Red" = "#F21A00"
  )
  
  globals$wes_fantasticfox1_colors <- c(
    "Orange" = "#DD8D29",
    "Yellow" = "#E2D200",
    "Teal" = "#46ACC8",
    "Burnt_Orange" = "#E58601",
    "Brown" = "#B40F20"
  )
  
  globals$wes_moonrise1_colors <- c(
    "Yellow" = "#F3DF6C",
    "Dark_Yellow" = "#CEAB07",
    "Light_Gray" = "#D5D5D3",
    "Dark_Gray" = "#24281A"
  )
  
# cluster colors and labels -------

  ## bladder cancer clusters 
  ## color-blind safe: 
  ## https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
  
  globals$cluster_levels <- c("#1", "#2", "#3")
  
  globals$cluster_colors <- c("#1" = "#009E73",
                              "#2" = "#0072B2",
                              "#3" = "#D55E00",
                              "unassigned" = "gray60")
  
  ## classification systems
  
  globals$system_colors <-
    c("clust_id" = "coral2",
      "NMIBC_class" = "aquamarine3",
      "consensusClass" = "steelblue3")
  
  globals$system_levels <-
    c("clust_id", "NMIBC_class", "consensusClass")
  
  globals$system_labels <-
    c("clust_id" = "bladder cancer clusters",
      "NMIBC_class" = "UROMOL NMIBC classes",
      "consensusClass" = "consensus MIBC classes")
  
  globals$system_shapes <-
    c("clust_id" = 21,
      "NMIBC_class" = 23,
      "consensusClass" = 25)
  
  ## consensus MIBC classes and consensus UROMOL classes 
  ## for NMIBC
  
  globals$consensus_colors <-
    c("LumP" = "coral4",
      "LumNS" = "coral2",
      "LumU" = "darkgoldenrod4",
      "Stroma-rich" = "chartreuse4",
      "Ba/Sq" = "dodgerblue3",
      "NE-like" = "plum4")
  
  globals$consensus_levels <- 
    c("LumP", "LumNS", "LumU", "Stroma-rich", "Ba/Sq", "NE-like")
  
  globals$uromol_colors <-
    c("1" = "aquamarine4", 
      "2a" = "steelblue2", 
      "2b" = "pink3", 
      "3" = "firebrick4")
  
  globals$uromol_levels <- c("1", "2a", "2b", "3")
  
  ## regulation colors and levels
  
  globals$regulation_colors <- 
    c("upregulated" = "firebrick", 
      "downregulated" = "steelblue", 
      "ns" = "gray60")
  
  globals$regulation_levels <- 
    c("upregulated", "downregulated", "ns")
  
# vertex attributes for similarity networks --------
  
  globals$attr_df <- 
    map2_dfr(globals[c("cluster_levels", "uromol_levels", "consensus_levels")], 
             globals$system_levels, 
             ~tibble(name = .x, 
                     system = .y)) %>% 
    mutate(system = factor(system, globals$system_levels))
  
# END ---------
