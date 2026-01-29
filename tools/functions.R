# Accessory functions for the project

# differential expression analyses ---------

  format_anova <- function(anova, 
                           annotation, 
                           p_cutoff = 0.05, 
                           etasq_cutoff = 0.06) {
    
    ## consistent formatting of one-way ANOVA results: 
    ## appending with Entrez ID and identificaiton of 
    ## significant effects with cutoffs of adjusted p and effect size
    
    anova %>% 
      left_join(annotation, by = "variable") %>% 
      mutate(regulation = ifelse(p_adjusted < p_cutoff & etasq >= etasq_cutoff, 
                                 "regulated", "ns"), 
             regulation = factor(regulation, c("regulated", "ns")))

  }
  
  format_posthoc <- function(test, 
                             annotation, 
                             anova_significant, 
                             p_cutoff = 0.05, 
                             regulation_levels = globals$regulation_levels) {
    
    ## consistent formatting of post-hoc one-sample T test for 
    ## effects in a subset as compared with the cohort average
    ## annotation with Entrez ID and identification of significant 
    ## effects with a text vector of variables found significant in one-way 
    ## ANOVA and a cutoff of adjusted p value
    
    test %>% 
      left_join(annotation, 
                by = "variable") %>% 
      mutate(anova_significant = ifelse(variable %in% anova_significant, 
                                        "yes", "no"), 
             anova_significant = factor(anova_significant, c("no", "yes"))) %>% 
      mutate(regulation = ifelse(anova_significant == "no" | p_adjusted >= 0.05, 
                                 "ns", 
                                 ifelse(deviation_center > 0, 
                                        "upregulated", 
                                        "downregulated")), 
             regulation = factor(regulation, regulation_levels))
    
  }
  
  format_roc <- function(roc, annotation, auc_cutoff = 0.714) {
    
    ## consistent formatting of ROC testing for markers of molecular 
    ## subsets
    
    roc %>% 
      mutate(marker = ifelse(auc >= auc_cutoff, "yes", "no"), 
             maker = factor(marker, c("yes", "no"))) %>% 
      arrange(-auc) %>% 
      left_join(annotation, by = "variable")
    
  }
  
  find_significant <- function(test_lst, split_fct) {
    
    ## extracts significant effects from post-hoc test results in a list
    
    test_lst %>% 
      map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
      map(function(x) if(nrow(x) == 0) NULL else x) %>% 
      compact %>% 
      map(blast, all_of(split_fct)) %>% 
      transpose %>% 
      map(compact) %>% 
      map(map, blast, regulation) %>% 
      map(transpose) %>% 
      map(map, map, ~.$variable)
    
  }
  
  find_markers <- function(roc_lst, split_fct) {
    
    ## extraction of markers of molecular subsets from 
    ## a list with ROC analyses
    
    roc_lst %>% 
      map(filter, marker == "yes") %>% 
      map(blast, all_of(split_fct)) %>% 
      transpose %>% 
      map(map, ~.x$variable)
    
  }
  
  count_significant <- function(test_lst, split_fct) {
    
    ## counts significant effects in a list of post-hoc results
    
    test_lst %>% 
      map(count, .data[[split_fct]], regulation) %>% 
      map(group_by, .data[[split_fct]]) %>% 
      map(mutate, 
          n_total = sum(n), 
          percent = n/n_total * 100) %>% 
      map(ungroup) %>% 
      map(filter, regulation %in% c("upregulated", "downregulated")) %>% 
      compress(names_to = "cohort")
    
  }
  
  count_cmm_significant <- function(x, 
                                    split_fct, 
                                    universe, 
                                    split_levels, 
                                    regulation_levels = globals$regulation_levels) {
    
    ## counts common significant features
    
    x %>% 
      map(map_dbl, length) %>% 
      map(compress, names_to = "regulation", values_to = "n") %>% 
      compress(names_to = split_fct) %>% 
      mutate(cohort = "common", 
             !!split_fct := factor(.data[[split_fct]], split_levels), 
             regulation = factor(regulation, regulation_levels), 
             n_total = length(universe), 
             percent = n/n_total * 100)
    
  }
  
  find_overlaps <- function(x, 
                            pairs, 
                            as_list = FALSE, 
                            regulation_levels = globals$regulation_levels) {
    
    ## extracts shared features from a list of character vectors (x), 
    ## the overlaps to check are defined by pairs: a list with vectors of 
    ## names of elements in x
    
    ## pair names and identifiers
    
    pair_names <- map(pairs, paste, collapse = "|")
    
    pairs <- set_names(pairs, pair_names)
    
    pair_df <- map(pairs, 
                   ~tibble(element1 = .x[[1]], 
                           element2 = .x[[2]]))
    
    pair_df <- compress(pair_df, names_to = "pair_id")
    
    ## overlaps
    
    res <- map(pairs, ~x[.x])
    
    res <- map(res, reduce, intersect)
    
    if(as_list) return(res)
    
    res <- map(res, 
               ~tibble(feature = .x))
    
    res <- compress(res, names_to = "pair_id")
    
    res <- left_join(res, pair_df, by = "pair_id")
    
    res <- mutate(res, 
                  variable = stri_split_fixed(feature, 
                                              pattern = "|", 
                                              simplify = TRUE)[, 1], 
                  regulation = stri_split_fixed(feature, 
                                                pattern = "|", 
                                                simplify = TRUE)[, 2], 
                  regulation = factor(regulation, regulation_levels))
    
    return(res)
    
  }
  
# Euler/Venn plots ---------
  
  plot_euler <- function(combinations, 
                         show_quantities = TRUE, 
                         show_labels = TRUE, 
                         alpha = 0.5, 
                         txt_size = 2.75, 
                         plot_title = NULL, 
                         plot_subtitle = NULL, 
                         cust_theme = globals$net_theme, ...) {
    
    ## a two-set Euler plot in a ggplot format
    ## tributes to https://gist.github.com/danlooo/
    
    ## plotting data ----------
    
    plot_lst <- eulerr::euler(combinations) %>% 
      plot(quantities = show_quantities) %>% 
      pluck("data")
    
    ellipse_tbl <- plot_lst$ellipses %>% 
      as_tibble(rownames = "Set")
    
    center_tbl <- plot_lst$centers %>% 
      mutate(label = NA)
    
    if(show_labels) {
      
      if(show_quantities) {
        
        center_tbl <- center_tbl %>% 
          mutate(label = ifelse(is.na(labels), 
                                quantities, 
                                paste(labels, quantities, sep = "\n")))
        
      } else {
        
        center_tbl <- center_tbl %>% 
          mutate(label = labels)
        
      }
      
    } else {
      
      if(show_quantities) {
        
        center_tbl <- center_tbl %>% 
          mutate(label = quantities)
        
      }
      
    }
    
    ## the plot -----------
    
    venn_plot <- tibble() %>% 
      ggplot() + 
      ggforce::geom_ellipse(data = ellipse_tbl, 
                            aes(x0 = h, 
                                y0 = k, 
                                a = a, 
                                b = b, 
                                angle = 0, 
                                fill = Set), 
                            alpha = alpha) + 
      cust_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle)
    
    if(show_quantities | show_labels) {
      
      venn_plot <- venn_plot + 
        geom_text(data = center_tbl, 
                  aes(x = x, y = y, label = label), 
                  size = txt_size)
      
    }
    
    venn_plot
    
  }
  
  euler_from_df <- function(x, 
                            subtract = TRUE, 
                            set1_n_var, 
                            set2_n_var, 
                            inter_n_var, 
                            set1_label, 
                            set2_label, 
                            suffix_label = NULL, 
                            plot_names = NULL, 
                            palette = globals[c("cluster_colors", 
                                                "consensus_colors", 
                                                "uromol_colors")] %>% 
                              reduce(c), 
                            fill_name = "molecular subset", 
                            ...) {
    
    ## makes a series of Euler plots with set counts and labels in a data frame
    ## subtract: counts of the intersection are subtracted from counts of set 1 
    ## and set 2
    
    ## plotting data -------
    
    if(subtract) {
      
      x[[set1_n_var]] <- x[[set1_n_var]] - x[[inter_n_var]]
      x[[set2_n_var]] <- x[[set2_n_var]] - x[[inter_n_var]]
      
    }
    
    combinations <- 1:nrow(x) %>% 
      map(function(idx) c(x[idx, set1_n_var][[1]], 
                          x[idx, set2_n_var][[1]], 
                          x[idx, inter_n_var][[1]]) %>% 
            set_names(as.character(x[idx, set1_label][[1]]), 
                      as.character(x[idx, set2_label][[1]]), 
                      paste(x[idx, set1_label][[1]], 
                            x[idx, set2_label][[1]], 
                            sep = "&")))
    
    if(!is.null(plot_names)) {
      
      plot_nm_vec <- as.character(x[[plot_names]])
      
      if(!is_null(suffix_label)) {
        
        plot_nm_vec <- paste(plot_nm_vec, 
                             x[[suffix_label]], 
                             sep = ".")
        
      }
      
      combinations <- set_names(combinations, plot_nm_vec)
      
    }
    
    ## plotting meta-data --------
    
    plot_titles <- paste(as.character(x[[set1_label]]), 
                         as.character(x[[set2_label]]), 
                         sep =  " and ")
    
    if(!is.null(suffix_label)) {
      
      plot_titles <- 
        paste(plot_titles, 
              as.character(x[[suffix_label]]), 
              sep = ", ")
      
    }
    
    ## plots -------
    
    list(combinations = combinations, 
         plot_title = plot_titles) %>% 
      pmap(plot_euler, ...) %>% 
      map(~.x + scale_fill_manual(values = palette, 
                                  name = fill_name))
    
  }
  
  euler_add_label <- function(x, 
                              label, 
                              x_offset = 1, 
                              txt_size = 2.75, 
                              txt_color = "black", 
                              hjust = 0, 
                              vjust = 0.5, 
                              ...) {
    
    ## adds text to an Euler plot: placed on the right side
    
    ## finding the text position
    
    set2_center_x <- x$layers[[1]]$data$h[[2]]
    set2_a <- x$layers[[1]]$data$a[[2]]
    
    txt_x <- set2_center_x + set2_a + x_offset
    
    x + 
      annotate("text", 
               label = label, 
               x = txt_x, 
               y = 0, 
               hjust = hjust, 
               vjust = vjust, 
               size = txt_size, 
               color = txt_color, ...)
    
  }
  
# General utilities --------
  
  mtx2long <- function(x, 
                       row_var = "variable1", 
                       col_var = "variable2", 
                       value_var = "value") {
    
    ## converts a matrix to long-format data frame: 
    ## variable 1: factors in rows, 
    ## variable 2: factors in columns
    
    x <- as.data.frame(x)
    
    cols_names <- colnames(x)
    
    x <- rownames_to_column(x, row_var)
    
    pivot_longer(x, 
                 cols = all_of(cols_names), 
                 names_to = col_var, 
                 values_to = value_var)
    
  }
  
  pad_missing <- function(x, variables, value = NA_real_) {
    
    missing_vars <- setdiff(variables, names(x))
    
    if(length(missing_vars) == 0) return(x)
    
    for(i in missing_vars) {
      
      x[[i]] <- value
      
    }
    
    return(x)
    
  }
  
# END --------
