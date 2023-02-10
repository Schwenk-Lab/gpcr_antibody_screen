### GPCR analysis miscellaneous functions
# Leo Dahl

# Validation of Abs
ab_validation <- function(long_dat, prot,
                          n_sd=12, excl_posctrl=T) {
  gpcrs <- long_dat %>%
    # Skip CALCRL.299 as it has no measurements for non-control CALCRL
    filter(gene_name_id != "CALCRL.299") %>%
    pull(gene_name_id) %>%
    unique() %>%
    sort()
  
  det_res <- sapply(gpcrs, function(ab) {
    d <- long_dat %>% filter(gene_name_id == ab, is.na(exclude))
    if (excl_posctrl) {d <- d %>% filter(gpcr_code != "positive IK19")}
    
    # Get density, peak and distance from peak position for computing sd and cutoff
    dens <- density(d$value)
    peak <- dens$x[which.max(dens$y)]

    pdist <- abs(d$value - peak)
    # Estimate expected proportion negative samples and use to get get values around peak to get variation of
    neg_prop <- 1 - (sum(d$target == "Target") / sum(d$target != "Target"))

    # SD or MAD depending on input argument
    # peak_var <- mad(d$value[pdist <= quantile(pdist, probs = neg_prop)])
    peak_var <- sd(d$value[pdist <= quantile(pdist, probs = neg_prop)])

    zsc_thr <- peak + n_sd * peak_var

    # Don't count positive control CALCRL as on-target as it's a different construct
    d[d$gpcr_code == "positive IK19", "target"] <- "Other"
    
    # GPCR expression detection
    gn <- prot %>% filter(gene_name_id == ab) %>% pull(gene_name) %>% unique()
    expr_detected <- !any("ns" %in% (long_dat %>% filter(target == "Target", gene_name == gn) %>%
                                       pull(gpcr_expr_ns_flag_capture_1d4_detect)))
    expr_detected <- ifelse(expr_detected, "supported", "uncertain")
    
    # Target detection
    target_detected <- any(d$target == "Target" & d$value >= zsc_thr)
    n_targets_detected <- sum(d$target == "Target" & d$value >= zsc_thr)
    n_targets_present <- sum(d$target == "Target")
    all_target_detected <- n_targets_detected == n_targets_present
    # Check if mean signal is above threshold
    mean_z <- d %>%
      filter(target == "Target") %>%
      pull(value) %>%
      mean()
    mean_z_pass <- mean_z > zsc_thr
    
    # Off-target detection
    # Special case for CALCRL. Since the control and target are the same GPCR, don't consider pos ctrl CALCRL binding to be cross-reactive for CALCRL binders
    # Reclassify the posctrl as "Target" now that the target detection is done
    if (grepl("CALCRL", ab)) {
      d[d$gpcr_code == "positive IK19", "target"] <- "Target"
    }
    
    other_detected <- any(d$target == "Other" & d$value >= zsc_thr)
    n_others_detected <- sum(d$target == "Other" & d$value >= zsc_thr)
    other_name_unique <- d %>%
      filter(target == "Other" & value >= zsc_thr) %>%
      pull(gpcr) %>%
      unique() %>%
      sort() %>%
      paste(., collapse = ", ")
    
    # Pass or fail
    pass <- all(c(mean_z_pass, !other_detected))
    
    df_out <- data.frame(
      "gene_name_id" = ab,
      "antibody" = unique(prot[prot$gene_name_id == ab, "antibody"]),
      "gpcr_expressed" = expr_detected,
      "signal_threshold" = signif(zsc_thr, 4),
      "mean_signal" = signif(mean_z, 4),
      "mean_signal_above_thr" = mean_z_pass,
      "pass" = pass,
      "target_detected" = target_detected,
      "n_targets_det" = n_targets_detected,
      "targets_total" = n_targets_present,
      "all_targets_det" = all_target_detected,
      "other_detected" = other_detected,
      "n_others_det" = n_others_detected,
      "other_names_unique" = other_name_unique
    )
    
    return(df_out)
  }, simplify = F, USE.NAMES = T) %>%
    rbindlist() %>%
    filter(!duplicated(gene_name_id))
  
}

# Making upset plot for Ab validation summary, colour just the pass bar
mk_upset <- function(x, sets,
                     main_colour = RColorBrewer::brewer.pal(3, "Set2")[1],
                     sub_colour = "grey23", matrix_lab_pos = "left",
                     tag = NULL, show_set_size = upset_set_size()) {
  # x, data frame to feed the upset function (containing 0 and 1). Should have a pass column named "Pass" (also with 0s and 1s)
  # sets, sets for the upset function
  # main_colour, colour to use for the pass bar
  # sub_colour, colour to use for the other bars
  
  # Colour bars based on pass/fail, include a legend
  x$Pass <- case_when(x$Pass == 1 ~ "Supported",
                      x$Pass == 0 ~ "Uncertain")
  plt <- ComplexUpset::upset(x, intersect = sets, name = "", set_sizes = show_set_size,
               base_annotations = list(
                 "Intersection size" = intersection_size(mapping = aes(fill = Pass), counts = T,
                                                         text = list(size = 6)) +
                   scale_fill_manual(values = c("Supported" = main_colour, "Uncertain" = sub_colour)) +
                   labs(fill = NULL, tag = tag) + 
                   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         axis.line = element_line(size = 1.5), axis.text = element_text(size = 20, colour = "black"),
                         axis.ticks = element_line(colour = "black", size = 1.5) ,axis.title = element_text(size = 20),
                         legend.text = element_text(size = 20), plot.tag = element_text(size = 24))),
               matrix = intersection_matrix() +
                 scale_y_discrete(position = matrix_lab_pos),
               themes = upset_modify_themes(
                 list("intersections_matrix" = theme(text = element_text(size = 24)))
               ))
  
  return(plt)
}

# Making a histogram for antibody validation
mk_hist <- function(x, binfo) {
  hist_df <- x %>%
    mutate(gpcr = binfo[match(gene_name_id, binfo$gene_name_id), "gene_name"]) %>%
    group_by(gpcr) %>%
    summarise(s = sum(pass, na.rm = T))
  
  ggplot(hist_df) +
    geom_histogram(aes(x = s), bins = length(0:max(hist_df$s))) +
    theme_classic(20) +
    labs(x = "Number of validated Abs", y = "Number of GPCRs") +
    scale_x_continuous(breaks = 0:max(hist_df$s)) +
    theme(axis.line = element_line(size = 1.5),
          axis.text = element_text(colour = "black", size = 20),
          axis.ticks = element_line(colour = "black", size = 1.5))
  
}

# Function for making GPCR expression plot
# Can swap p-value for asterisks instead
mk_gpcr_expr_plot <- function(dat_in, stars = F, log_in = F) {
  gpcr_expr <- dat_in$long_dat_all %>%
    filter(sample_type != "empty",
           !grepl("p11", sba, ignore.case = T),
           gene_name_id != "CALCRL.299") %>%
    filter((grepl("FZD", gpcr) & gene_name == "M anti-HA") |
             (!grepl("FZD", gpcr) & gene_name == "M anti-FLAG") |
             (sample_type == "mock" & gene_name == "M anti-HA")) %>%
    mutate(gexpr = case_when(grepl("mock", gpcr, ignore.case = T) ~ "Mock",
                             gpcr_expr_ns_flag_capture_1d4_detect == "ns" ~ "n.s.",
                             sample_type == "mock" ~ "Mock",
                             T ~ "Expressed")) 
  
  plt_out <- gpcr_expr %>%
    group_by(gene_name) %>%
    nest() %>%
    mutate(plt = map(data, ~ {
      plt_dat <- .x %>%
        select(subfamily, value) %>%
        # Rename subfamilies to not have NA and to bunch together the GSA group
        mutate(subfamily = case_when(is.na(subfamily) ~ "mock",
                                     subfamily %in% c("adhesion", "glutamate", "secretin") ~ "GSA",
                                     T ~ subfamily),
               subfamily = as.factor(subfamily)) # Factor to work with dunnett test function
      
      # ANOVA
      anova_p <- summary(aov(value ~ subfamily, plt_dat))[[1]]["subfamily", "Pr(>F)"]
      
      # Dunnett's test for multiple comparisons if the ANOVA p < 0.05
      if (anova_p < 0.05) {
        dunnett_p <- DunnettTest(x = plt_dat$value, g = plt_dat$subfamily, control = "mock")
        dunnett_p <- as.data.frame(dunnett_p$mock) %>%
          rownames_to_column("subfamily") %>%
          select(subfamily, pval) %>%
          # Remove control group name, rename to match names on plot axis
          mutate(subfamily = str_remove(subfamily, "-.*"))
        
        # The test gives 0 for p-values lower than 2e-16, if not asterisks swap zeroes (what is given for values < 2e-16) with "2e-16"
        if (!stars) {
          dunnett_p$pval <- case_when(dunnett_p$pval < 2e-16 ~ "<2e-16",
                                      T ~ as.character(signif(dunnett_p$pval, 3)))
        } else {
          # Substitute numbers for asterisks with borders at 0.0001, 0.001, 0.01 and 0.05
          dunnett_p$pval <- symnum(dunnett_p$pval, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                   symbols = c("****", "***", "**", "*", "n.s."), corr = F, legend = F) %>% as.character()
        }
      }
      
      # Adjust subfamily labels
      # Change from factor to character vector, otherwise case_when() gives an error
      plt_dat$subfamily <- as.character(plt_dat$subfamily)
      
      # Log values if not already log in-values
      if (!log_in) {plt_dat$value <- log2(plt_dat$value)}
      
      # Reorder x-axis labels
      subfam_order <- c("mock", "alpha", "beta", "gamma", "delta", "GSA", "other", "frizzled")
      plt_dat$subfamily <- factor(plt_dat$subfamily, levels = subfam_order)
      
      plt <- ggplot(plt_dat, aes(x = subfamily, y = value)) +
        geom_beeswarm(colour = "grey70") +
        # stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
        #              size = 0.4, width = ifelse(gene_name == "M anti-FLAG", 0.5, 0.15), colour = "black") +
        stat_summary(fun = median, fun.min = median, fun.max = median,
                     geom = "errorbar", colour = "grey40", width = ifelse(gene_name == "M anti-FLAG", 0.2, 0.15), size = 1.2) +
        stat_summary(fun = median,
                     fun.min = function(a) quantile(a, 0.25),
                     fun.max = function(a) quantile(a, 0.75),
                     geom = "errorbar", colour = "grey40", width = ifelse(gene_name == "M anti-FLAG", 0.15, 0.1), size = 0.9) +
        labs(tag = ifelse(gene_name == "M anti-FLAG", "A", "B"), x = "Subfamily", y = paste0(ifelse(gene_name == "M anti-FLAG", "FLAG", "HA"), "-capture\nlog2(MFI) [AU]")) +
        scale_x_discrete(labels = c("mock" = "Mock",
                                    "alpha" = expression(alpha),
                                    "beta" = expression(beta),
                                    "gamma" = expression(gamma),
                                    "delta" = expression(delta),
                                    "GSA" = "GSA",
                                    "other" = "Other",
                                    "frizzled" = "Frizzled")) +
        ylim(4, 15) +
        theme_classic(14) +
        theme(axis.text = element_text(size = 14, colour = "black"),
              axis.ticks = element_line(colour = "black", size = 1),
              axis.line = element_line(size = 1))
      
      # Add p-values if relevant
      if (anova_p < 0.05) {
        # Get y range from plot to place p-values
        y_range <- ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range
        y_pos <- max(y_range) + 0.05 * diff(y_range) # Place slightly above max y
        
        plt <- plt +
          geom_text(data = dunnett_p, aes(label = pval), y = y_pos, size = 6) +
          coord_cartesian(clip = "off")
      }
      
      return(plt)
    })) %>%
    pull(plt)
  
  return(plt_out)
}


# Look at GPCR expression per GPCR (not per subfamily)
per_gpcr_expr <- function(dat_in) {
  gpcr_expr <- dat_in$long_dat_all %>%
    filter(sample_type != "empty",
           !grepl("p11", sba, ignore.case = T),
           gene_name_id != "CALCRL.299",
           gpcr_code != "positive IK19") %>%
    filter((grepl("FZD", gpcr) & gene_name == "M anti-HA") |
             (!grepl("FZD", gpcr) & gene_name == "M anti-FLAG") |
             (sample_type == "mock" & gene_name == "M anti-HA")) %>%
    mutate(gexpr = case_when(grepl("mock", gpcr, ignore.case = T) ~ "Mock",
                             gpcr_expr_ns_flag_capture_1d4_detect == "ns" ~ "n.s.",
                             sample_type == "mock" ~ "Mock",
                             T ~ "Expressed"))
  
  # Perform ANOVA + Dunnett test separately for each capture
  gpcr_expr %>%
    group_by(gene_name) %>%
    nest() %>%
    mutate(p_anova = map(data, ~ {
      summary(aov(value ~ gpcr, .x))
    }))
  
}

# Plot comparing passing Abs to Abs with offtarget binding for different numbers of sd
pass_vs_sd <- function(dat_list, n_sd) {
  plt_dat <- sapply(n_sd, function(x) {
    d <- ab_validation(dat_list$long_dat, dat_list$prot, x, T)
    
    n_pass <- sum(d$pass, na.rm = T)
    n_xrxn <- sum(d$other_detected, na.rm = T)
    
    df_out <- data.frame("n_sd" = x,
                         "On_target" = n_pass,
                         "Crossreactive" = n_xrxn)
  }, simplify = F) %>% rbindlist()
  
  plt_dat %>%
    pivot_longer(cols = -n_sd) %>%
    ggplot() +
    geom_line(aes(x = n_sd, y = value, colour = name), size = 3) +
    scale_colour_brewer(palette = "Set2", direction = -1) +
    labs(x = "Number of SD", y = "Number of Ab", colour = "Type") +
    theme_classic() +
    theme(axis.line = element_line(size = 1.5),
          axis.text = element_text(size = 20, colour = "black"),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          axis.ticks = element_line(size = 1.5, colour = "black"))
}

