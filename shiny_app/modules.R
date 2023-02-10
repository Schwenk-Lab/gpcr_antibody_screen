### --- GPCR Shiny app Modules --- ###

## List of GPCRs for selection in analysis
gpcr_select_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("gpcr_select")),
    uiOutput(ns("uniprot_select"))
  )
}

gpcr_select_server <- function(id, data_in) {
  # data_in, dataset containing the shorter long data structure
  moduleServer(id, function(input, output, session) {
    
    # GPCR names to choose from
    gpcrs <- data_in$abval_dat %>%
      distinct(gpcr) %>%
      filter(!gpcr %in% c("empty", "mock")) %>%
      pull(gpcr) %>%
      sort() %>%
      `names<-`(., NULL)
    # UniProt IDs in same order as GPCRs
    uniprots <- data_in$abval_dat[match(gpcrs, data_in$abval_dat$gpcr), "uniprot", drop = T]
    # Matching GPCR and uniprot for syncing
    gpcr_uniprot <- gpcrs %>% `names<-`(., uniprots)
    
    # Make GPCR input UI with GPCRs from data
    output$gpcr_select <- renderUI({
      selectInput(NS(id, "choose_gpcr"), "Select GPCR", gpcrs, selected = "ADORA1")
    })
    
    output$uniprot_select <- renderUI({
      selectInput(NS(id, "choose_uniprot"), "Select UniProt ID", uniprots, selected = "P30542")
    })
    
    # Sync the inputs
    observeEvent(input$choose_uniprot, {
      updateSelectInput(session, "choose_gpcr", "Select GPCR", gpcrs, gpcr_uniprot[input$choose_uniprot])
    })
    observeEvent(input$choose_gpcr, {
      updateSelectInput(session, "choose_uniprot", "Select UniProt ID", uniprots, names(gpcr_uniprot)[gpcr_uniprot == input$choose_gpcr])
    })
    
    return(reactive(input$choose_gpcr))
  })
}

### Antibody validation ###
## Beeswarm plot for antibody validation
abval_beeswarm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Antibody validation beeswarms (1D4 detection)"),
    h4("Plot customisation"),
    # Plot customisation for Ab validation beeswarm plots, allow resizing of points
    sliderInput(ns("point_size"), "Point size", min = 1, max = 10, value = 2, step = 0.5),
    hr(),
    plotOutput(ns("beeswarm"), height = 600)
  )
}

abval_beeswarm_server <- function(id, data_in, gpcr) {
  # data_in: Data set containing the abval_dat data frame and the Ab information
  # gpcr: Reactive selected GPCR from the gpcr_select module
  moduleServer(id, function(input, output, session) {
    output$beeswarm <- renderPlot({
      # Prevent errors from non-existent GPCR selection
      gpcr_in <- req(gpcr())
      
      # To deal with multitarget Abs, get the Abs targeting the selected GPCR beforehand
      abs <- data_in$prot %>%
        # Check each Ab if its targets contain the selected GPCR
        filter(map_lgl(data_in$prot$gene_name, function(x) {
          gpcr_in %in% (strsplit(x, ";") %>% unlist())
        })) %>%
        pull(gene_name_id)
      
      # Keep only selected GPCR
      d <- data_in$abval_dat %>%
        # Regex to deal with multitarget Abs
        filter(gene_name_id %in% abs) %>%
        # Exclude failed samples
        filter(is.na(exclude))
      
      if (nrow(d) == 0) {
        return(ggplot() + labs(title = paste0("No antibodies for ", gpcr_in)))
      } else {
        p <- sapply(sort(unique(d$gene_id_ab)), function(ab) {
          d_temp <- d %>% filter(gene_id_ab == ab)
          ab_name <- str_split(ab, "\\.") %>% unlist()
          ab_name <- paste0(ab_name[1], "_", ab_name[3])
          
          gni <- unique(d_temp$gene_name_id)
          thr <- data_in$prot %>% filter(gene_name_id == gni) %>% pull(threshold)
          
          ggplot(d_temp, aes(x = 0, y = value_rz, colour = target, fill = target)) +
            geom_beeswarm(size = input$point_size, shape = 21, alpha = 0.6, priority = "none") +
            geom_hline(yintercept = thr, linetype = "dashed") +
            scale_colour_manual(values = c("Target" = "blue", "Other" = "grey")) +
            scale_fill_manual(values = c("Target" = "blue", "Other" = "grey")) +
            scale_x_continuous(labels = NULL, breaks = NULL) +
            labs(x = NULL, y = "Robust Z-score",
                 title = ab_name, colour = NULL, fill = NULL) +
            theme_bw() +
            theme(axis.line = element_line(size = 1),
                  plot.title = element_text(size = 15),
                  axis.text = element_text(size = 12),
                  legend.text = element_text(size = 12))
          
        }, simplify = F)
        
        return(wrap_plots(p, guides = "collect"))
      }
    })
  })
}


## Table of antibody information
ab_info_ui <- function(id) {
  ns <- NS(id)
  DTOutput(ns("ab_tbl"))
}

ab_info_server <- function(id, data_in, gpcr) {
  # data_in, data set containing antibody information
  # gpcr, reactive string of the selected GPCR name, or non-reactive input to show all
  moduleServer(id, function(input, output, session) {
    output$ab_tbl <- renderDT({
      
      # Show all GPCRs or just one selected one
      if (!is.reactive(req(gpcr))) {
        # Remove Abs with non-GPCR targets
        dt_in <- data_in$prot %>%
          filter(!gene_name %in% c("M anti-FLAG", "M anti-HA")) %>%
          filter(!is.na(crossreact_names_unique)) %>%
          arrange(gene_name)
      } else {
        gpcr_in <- req(gpcr())
        # Get info for Abs targeting selected GPCR
        dt_in <- data_in$prot %>%
          # Check each Ab to see if its targets contain the selected GPCR
          filter(map_lgl(data_in$prot$gene_name, function(x) {
            gpcr_in %in% (strsplit(x, ";") %>% unlist())
          }))
      }
      
      dt_in %<>%
        filter(gene_name_id != "CALCRL.299") %>%
        select(antibody, gene_name, rrid, pass, crossreact, crossreact_names_unique) %>%
        # Make annotation consistent across columns
        mutate(rrid = case_when(rrid == "" ~ "None", T ~ rrid),
               pass = case_when(pass ~ "Yes", T ~ "No"),
               crossreact = case_when(crossreact ~ "Yes", T ~ "No"),
               crossreact_names_unique = case_when(crossreact_names_unique == "" ~ "None",
                                                   T ~ crossreact_names_unique)) %>%
        rename("GPCR" = "gene_name", "Antibody" = "antibody", "RRID" = "rrid",
               "On-target" = "pass","Crossreacts" = "crossreact",
               "Other target" = "crossreact_names_unique") %>%
        distinct()
      
      datatable(dt_in, rownames = F, options = list(dom = ifelse(!is.reactive(req(gpcr)), "ltiprf", "tr"),
                                                    pageLength = ifelse(!is.reactive(req(gpcr)), 10, nrow(dt_in))))
    })
  }) 
}

## Expression plot for selected target
expr_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("GPCR expression"),
    plotOutput(ns("plt"), width = "90%")
  )
}

expr_plot_server <- function(id, data_in, gpcr) {
  # data_in: Data containing the expr_dat long data (for FLAG capture)
  # gpcr: reactive selected GPCR for Ab validation
  moduleServer(id, function(input, output, session) {
    output$plt <- renderPlot({
      # Prevent error from non-existent GPCR selection
      gpcr_in <- req(gpcr())
      
      d <- data_in$expr_dat %>%
        # FLAG capture, 1D4 detection
        # If FZD, capture with HA and look at only FZD, otherwise capture with FLAG and look at all except FZD (due to flipped FLAG/HA tags)
        filter(if (grepl("FZD", gpcr_in)) {gene_name == "M anti-HA" & (grepl("FZD", gpcr) | gpcr %in% c("empty", "mock"))}
               else {gene_name == "M anti-FLAG" & !grepl("FZD", gpcr)}) %>%
        filter(gpcr_code != "positive IK19") %>%
        mutate(value = log2(value_r)) %>%
        # Make separate densities for mock + empty and GPCRs
        mutate(grp = case_when(gpcr %in% c("mock", "empty") ~ "Negative control",
                               T ~ "GPCR"),
               grp = factor(grp, levels = c("Negative control", "GPCR")))
      
      # Whether the GPCR was deemed significantly expressed
      gpcr_expr <- !("ns" %in% (d %>% filter(gpcr == gpcr_in) %>% pull(gpcr_expr_ns_flag_capture_1d4_detect)))
      expr_col <- ifelse(gpcr_expr, brewer.pal(3, "Set2")[1], brewer.pal(3, "Set2")[2])
      
      ggplot(d %>% select(value, grp)) +
        geom_density(aes(x = value, y = ..scaled.., colour = grp), size = 1.5, adjust = 0.5) +
        geom_vline(data = filter(d, gpcr == gpcr_in), aes(xintercept = value, size = ifelse(gpcr_expr, "Supported", "Uncertain")),
                   linetype = "dotted", colour = "#4b99fe") +
        # Force in a shared legend for the density lines and vertical lines
        geom_rect(data = d %>% select(value, grp),
                  aes(alpha = factor(c("Negative control", "GPCR", rep("Target sample", nrow(d) - 2)),
                                     levels = c("Negative control", "GPCR", "Target sample")),
                      xmin = median(d$value), xmax = median(d$value), ymin = 0.5, ymax = 0.5)) +
        scale_alpha_manual("Sample type", values = 1:3,
                           guide = guide_legend(override.aes = list(size = 1.2,
                                                                    fill = "white",
                                                                    linetype = c(1, 1, 3),
                                                                    colour = c("#FC8D62", "#66C2A5", "#4b99fe")))) +
        guides(colour = "none") +
        # Force in a legend saying whether the expression is supported or not
        scale_size_manual(paste0(gpcr_in, "\nExpression"), values = 1,
                          guide = guide_legend(override.aes = list(colour = "white"), order = 1, label.position = "left",
                                               label.theme = element_text(size = 18, colour = expr_col))) +
        coord_cartesian(clip = "off") +
        scale_colour_brewer(palette = "Set2", direction = -1) +
        labs(x = "GPCR expression level", y = "Density", colour = "Sample type") +
        theme_classic(18) +
        theme(axis.line = element_line(size = 1))
      
    })
  })
}


## Table with GPCR information
gpcr_info_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("gpcr_tbl"))
}

gpcr_info_server <- function(id, data_in, gpcr) {
  # data_in, data set containing the long data
  # gpcr, reactive string of the selected GPCR
  moduleServer(id, function(input, output, session) {
    
    # Get one row of info for current GPCR
    current_gpcr <- reactive({
      data_in$abval_dat %>%
        filter(gpcr == req(gpcr())) %>%
        slice(1)
    })
    
    # GPCR table containing links to HPA, Ensembl, UniProt, different identifiers, info about family/subfamily
    output$gpcr_tbl <- renderUI({
      HTML(paste0(
        "<h1>", req(gpcr()), "</h1>",
        "<p>
        <a href='", paste0("https://www.proteinatlas.org/",
                           current_gpcr()$ensg, "-", req(gpcr())),
        "' class='gpcrlink' target='_blank' rel='noopener noreferrer'>HPA portal</a>
        <a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
                           current_gpcr()$ensg),
        "' class='gpcrlink' target='_blank' rel='noopener noreferrer'>Ensembl</a>
        <a href='", paste0("https://www.uniprot.org/uniprot/",
                           str_remove(current_gpcr()$uniprot, "[[:space:]]")), # Remove space at the end of uniprot IDs (if present), breaks link
        "' class='gpcrlink' target='_blank' rel='noopener noreferrer'>UniProt</a>
        </p>",
        "<table style='width:100%' class='gpcrinfo'>
          <tr>
            <td>Ensembl gene ID</td>
            <td>", current_gpcr()$ensg,"</td>
          </tr>
          <tr>
            <td>UniProt ID</td>
            <td>", current_gpcr()$uniprot, "</td>
          </tr>
          <tr>
            <td>Class</td>
            <td>", current_gpcr()$class, "</td>
          </tr>
          <tr>
            <td>Subfamily</td>
            <td>", current_gpcr()$subfamily, "</td>
          </tr>
          <tr>
            <td>Ligand Family</td>
            <td>", current_gpcr()$family, "</td>
          </tr>
        </table>"
      ))
    })
  })
}

### UpSet plot ###
## UpSet plot from antibody validation
abval_upset_ui <- function(id) {
  ns <- NS(id)
  plotOutput(ns("upset_plot"), height = "600px")
}

abval_upset_server <- function(id, data_in, sets) {
  moduleServer(id, function(input, output, session) {
    # data_in, data set containing protein table
    # sets, vector of selected sets to display
    
    # Make columns for upset plot
    d <- data_in$prot %>%
      filter(!gene_name %in% c("M anti-FLAG", "M anti-HA")) %>%
      filter(!duplicated(gene_name_id)) %>%
      select(gene_name_id, pass, detect_target, crossreact, gpcr_expr) %>%
      # Names to match selections
      rename("Pass" = "pass",
             "On-target detection" = "detect_target",
             "Crossreactivity" = "crossreact",
             "GPCR expression" = "gpcr_expr") %>%
      mutate(Pass = case_when(!is.na(`On-target detection`) ~ Pass),
             `GPCR expression` = case_when(`GPCR expression` == "supported" ~ T,
                                           `GPCR expression` == "uncertain" ~ F)) %>%
      distinct() %>%
      # Numeric columns instead of logical for the upset function to work
      mutate(across(where(is.logical), as.numeric)) %>%
      # Better names for pass/fail for legend
      mutate(Pass = case_when(Pass == 1 ~ "Supported", Pass == 0 ~ "Uncertain"))
    
    output$upset_plot <- renderPlot({
      
      ComplexUpset::upset(d %>% na.omit(), intersect = sets, name = "", set_sizes = F,
                          base_annotations = list(
                            "Intersection size" = intersection_size(text = list(size = 8), mapping = aes(fill = Pass), counts = T) +
                              ylab("Number of Antibodies") +
                              scale_fill_manual(values = c("Supported" = brewer.pal(3, "Set2")[1], "Uncertain" = "black")) + labs(fill = NULL) +
                              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    axis.line = element_line(size = 0.5), axis.text = element_text(size = 16),
                                    axis.title = element_text(size = 16), legend.text = element_text(size = 16))),
                          themes = upset_modify_themes(list(
                            "intersections_matrix" = theme(text = element_text(size = 18)),
                            "overall_sizes" = theme(axis.text.x = element_text(size = 14),
                                                    axis.title.x = element_text(size = 12))
                          )))
      
    })
  })
}

### Median signal heatmaps ###
## Median signal heatmaps, family and subfamily selection
fam_subfam_select_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("med_hm_fam_subfam"), "Show subfamilies or ligand families",
                 c("Subfamily", "Ligand family"), "Subfamily"),
    uiOutput(ns("fam_select"))
  )
}


fam_subfam_select_server <- function(id, data_in) {
  # data_in, data containing long format data with subfamily/family classifications
  
  moduleServer(id, function(input, output, session) {
    
    output$fam_select <- renderUI({
      if (input$med_hm_fam_subfam == "Subfamily") {
        options_out <- c("Rhodopsin (\u03B1)", "Rhodopsin (\u03B2)", "Rhodopsin (\u03B3)", "Rhodopsin (\u03B4)",
                         "GSAF", "Other + other SBA", "Other + all SBAs")
      } else if (input$med_hm_fam_subfam == "Ligand family") {
        options_out <- data_in$abval_dat %>%
          pull(family) %>%
          str_to_title() %>%
          unique() %>%
          sort() %>%
          na.omit()
        
      }
      
      selectInput(NS(id, "fam_select"),
                  ifelse(input$med_hm_fam_subfam == "Subfamily", "Select subfamily", "Select ligand family"),
                  options_out)
    })
    
    return(list(reactive(input$med_hm_fam_subfam), reactive(input$fam_select)))
  })
}

## Median signal heatmaps, customising the heatmap
# Relative heights controls for heatmaps with subplots, gives error when switching between heatmaps that have different numbers
# of subplots. But even with controls it does not look much better, and it becomes even more tricky to customise the plots. 
# The colourbars are also not scaling with the subplots, making it hard to know which bar is for which plot.
# The code is kept here in case the feature is requested
hm_custom_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    h4("Heatmap layout"),
    h5("Heatmap size (px)"),
    splitLayout(
      numericInput(ns("hm_height"), "Height", value = 800, step = 10),
      numericInput(ns("hm_width"), "Width", value = 800, step = 10)
    ),
    h5("Axis label sizes (0 to remove labels)"),
    splitLayout(
      numericInput(ns("ab_size"), "Antibody", value = 7, step = 1, min = 0),
      numericInput(ns("gpcr_size"), "GPCR", value = 10, step = 1, min = 0)
    ),
    checkboxInput(ns("transpose"), "Click to swap axes", value = F)
  )
}

hm_custom_server <- function(id) {
  # fam_selection, reactive, character of length 1 indicating if subfamily or family should be plotted. For conditional UI element
  moduleServer(id, function(input, output, session) {
    return(reactive(list(
      "height" = input$hm_height,
      "width" = input$hm_width,
      "ab_size" = input$ab_size,
      "gpcr_size" = input$gpcr_size,
      "transpose" = input$transpose
    )))
  })
}

## Median signal heatmaps, making the heatmap
med_hm_ui <- function(id) {
  ns <- NS(id)
  plotlyOutput(ns("heatmap"))
}

med_hm_server <- function(id, subfam_or_fam, fam_selection, data_in, param) {
  # subfam_or_fam, reactive, character of length 1 indicating if subfamily or family should be plotted
  # fam_selection, reactive, character of length 1 with the family or subfamily name, from fam_subfam_select module
  # data_in, data set
  # param, reactive list with heatmap layout parameters
  moduleServer(id, function(input, output, session) {
    output$heatmap <- renderPlotly({
      # If statement to prevent error when switching between subfamily/family. Since they both use the same fam_selection() reactive
      # variable, an error occurs when switching as the value of fam_selection() that already exists is used. So it will look for a 
      # subfamily within the ligand families or vice versa, not find it, and throw an error. Prevent by requiring fam_selection()
      # to be a subfamily when subfamily is selected or a ligand family when ligand family is selected
      if ((
        req(subfam_or_fam()) == "Subfamily" &
        req(fam_selection()) %in% c("Rhodopsin (\u03B1)", "Rhodopsin (\u03B2)", "Rhodopsin (\u03B3)", "Rhodopsin (\u03B4)",
                                    "GSAF", "Other + other SBA", "Other + all SBAs")
      ) | (
        req(subfam_or_fam()) == "Ligand family" &
        !req(fam_selection()) %in% c("Rhodopsin (\u03B1)", "Rhodopsin (\u03B2)", "Rhodopsin (\u03B3)", "Rhodopsin (\u03B4)",
                                     "GSAF", "Other + other SBA", "Other + all SBAs")
      )) {
        if (subfam_or_fam() == "Subfamily") {
          
          # Get data for plotting
          plt_dat <- data_in$abval_dat %>%
            mutate(sba = data_in$prot[match(gene_name_id, data_in$prot$gene_name_id), "sba"]) %>%
            filter(if (fam_selection() == "Other + other SBA") {sba == "other" & subfamily == "Other"}
                   else if (fam_selection() == "Other + all SBAs") {subfamily == "Other"}
                   else {subfamily == fam_selection()})
          
          # Title of plot
          plt_title <- paste0("Subfamily: ", fam_selection())
        } else if (subfam_or_fam() == "Ligand family") {
          plt_dat <- data_in$abval_dat %>%
            filter(family == fam_selection())
          
          plt_title <- paste0("Ligand family: ", fam_selection())
        }
        
        plt_dat %<>%
          filter(is.na(exclude)) %>%                          # Exclude samples marked for exclusion
          group_by(gene_name_id, gpcr) %>%
          mutate(median_signal = median(value_z)) %>%         # Calculate median signal per GPCR-Ab combination
          ungroup() %>%
          distinct() %>%
          # Axis labels containing target name and antibody name, sorting puts target and antibody on diagonal
          mutate(gpcr_ab = map_chr(gene_id_ab, function(x) {
            spl_name <- str_split(x, "\\.") %>% unlist()
            paste0(spl_name[1], "_", spl_name[3])
          })) %>%
          select(gpcr, gpcr_ab, median_signal)
        
        # Some ligand families have large gaps due to the GPCRs being assayed with different SBA panels
        # For such cases, make multiple heatmaps to skip the gaps
        if (fam_selection() %in% c("Complement Peptide", "Leukotriene", "Lysophospholipid", "Orphan", "Somatostatin")) {
          
          plt_out <- plt_dat %>%
            # Fill in all combinations to count missing values
            complete(gpcr, gpcr_ab) %>%
            # Count missing values per binder to make a plot per unique number of missing values
            group_by(gpcr_ab) %>%
            mutate(nmiss = sum(is.na(median_signal))) %>%
            ungroup() %>%
            # Remove missing values to not get in plot
            filter(!is.na(median_signal)) %>%
            nest_by(nmiss) %>%
            # Make plots
            mutate(plt = list(
              plot_ly(data = data,
                      x = ~ifelse(rep(param()$transpose, nrow(data)), gpcr_ab, gpcr),
                      y = ~ifelse(rep(param()$transpose, nrow(data)), gpcr, gpcr_ab),
                      z = ~median_signal, type = "heatmap", colors = viridis(100), xgap = 0.3, ygap = 0.3,
                      hoverinfo = "text",
                      hovertext = paste(
                        ifelse(rep(param()$transpose, nrow(data)), paste("Antibody:", data$gpcr_ab), paste("GPCR:", data$gpcr)),
                        ifelse(rep(param()$transpose, nrow(data)), paste("<br>GPCR:", data$gpcr), paste("<br>Antibody:", data$gpcr_ab)),
                        "<br>Median Z-score:", signif(data$median_signal, 4)
                      ),
                      height = param()$height, width = param()$width) %>%
                layout(title = plt_title,
                       xaxis = list(title = ifelse(param()$transpose, "Antibody", "GPCR"),
                                    nticks = ifelse(param()$transpose, length(unique(data$gpcr_ab)), length(unique(data$gpcr))),
                                    tickfont = list(size = ifelse(param()$transpose, param()$ab_size, param()$gpcr_size)),
                                    ticks = ifelse(param()$transpose,
                                                   ifelse(param()$ab_size <= 0, "", "outside"),
                                                   ifelse(param()$gpcr_size <= 0, "", "outside")),
                                    showticklabels = ifelse(param()$transpose,
                                                            ifelse(param()$ab_size <= 0, F, T),
                                                            ifelse(param()$gpcr_size <= 0, F, T))),
                       yaxis = list(title = ifelse(param()$transpose, "GPCR", "Antibody"),
                                    nticks = ifelse(param()$transpose, length(unique(data$gpcr)), length(unique(data$gpcr_ab))),
                                    tickfont = list(size = ifelse(param()$transpose, param()$gpcr_size, param()$ab_size)),
                                    ticks = ifelse(param()$transpose,
                                                   ifelse(param()$gpcr_size <= 0, "", "outside"),
                                                   ifelse(param()$ab_size <= 0, "", "outside")),
                                    showticklabels = ifelse(param()$transpose,
                                                            ifelse(param()$gpcr_size <= 0, F, T),
                                                            ifelse(param()$ab_size <= 0, F, T)))) %>%
                colorbar(title = "Median Z-score")
            )) %>% pull(plt) %>%
            subplot(nrows = length(.))
          
          return(plt_out)
        }
        
        # Text for mouse hover
        hvr_txt <- paste(ifelse(rep(param()$transpose, nrow(plt_dat)), paste("Antibody:", plt_dat$gpcr_ab), paste("GPCR:", plt_dat$gpcr)),
                         ifelse(rep(param()$transpose, nrow(plt_dat)), paste("<br>GPCR:", plt_dat$gpcr), paste("<br>Antibody:", plt_dat$gpcr_ab)),
                         "<br>Median Z-score:", signif(plt_dat$median_signal, 4))
        
        # The plot, lots of ifelse statements (and a few nested) because of the option to transpose the heatmap
        plt_out <- plot_ly(data = plt_dat,
                           x = ~ifelse(rep(param()$transpose, nrow(plt_dat)), gpcr_ab, gpcr),
                           y = ~ifelse(rep(param()$transpose, nrow(plt_dat)), gpcr, gpcr_ab),
                           z = ~median_signal, type = "heatmap", colors = viridis(100), xgap = 0.3, ygap = 0.3,
                           hoverinfo = "text", hovertext = hvr_txt,
                           height = param()$height, width = param()$width) %>%
          layout(title = plt_title,
                 xaxis = list(title = ifelse(param()$transpose, "Antibody", "GPCR"),
                              nticks = ifelse(param()$transpose, length(unique(plt_dat$gpcr_ab)), length(unique(plt_dat$gpcr))),
                              tickfont = list(size = ifelse(param()$transpose, param()$ab_size, param()$gpcr_size)),
                              ticks = ifelse(param()$transpose,
                                             ifelse(param()$ab_size <= 0, "", "outside"),
                                             ifelse(param()$gpcr_size <= 0, "", "outside")),
                              showticklabels = ifelse(param()$transpose,
                                                      ifelse(param()$ab_size <= 0, F, T),
                                                      ifelse(param()$gpcr_size <= 0, F, T))),
                 yaxis = list(title = ifelse(param()$transpose, "GPCR", "Antibody"),
                              nticks = ifelse(param()$transpose, length(unique(plt_dat$gpcr)), length(unique(plt_dat$gpcr_ab))),
                              tickfont = list(size = ifelse(param()$transpose, param()$gpcr_size, param()$ab_size)),
                              ticks = ifelse(param()$transpose,
                                             ifelse(param()$gpcr_size <= 0, "", "outside"),
                                             ifelse(param()$ab_size <= 0, "", "outside")),
                              showticklabels = ifelse(param()$transpose,
                                                      ifelse(param()$gpcr_size <= 0, F, T),
                                                      ifelse(param()$ab_size <= 0, F, T)))) %>%
          colorbar(title = "Median Z-score")
        
        return(plt_out)
      }
    })
  })
}


### Sibling antibody correlation ###
## Selecting antibody target (as GPCR name or UniProt ID)
target_select_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("select_target")),
    uiOutput(ns("select_uniprot"))
  )
}

target_select_server <- function(id, data_in) {
  # data_in, data containing long format data
  moduleServer(id, function(input, output, session) {
    
    # Ab targets (that have multiple Abs against them) and their corresponding uniprot ids to choose from
    ab_targets <- data_in$abval_dat %>%
      group_by(gene_name) %>%
      summarise(n_ab = length(unique(gene_name_id)), .groups = "keep") %>%
      filter(n_ab > 1) %>%
      pull(gene_name) %>%
      unique() %>%
      sort()
    uniprots <- data_in$abval_dat[match(ab_targets, data_in$abval_dat$gpcr), "uniprot", drop = T]
    # The targets with multiple GPCRs separated by ; get NA as ID, deal with separately
    # Their uniprot IDs are hard-coded here as some, like NPY4R2, are not included in any samples and thus have no listed uniprot IDs in the data
    uniprots[ab_targets == "HCAR2;HCAR3"] <- "Q8TDS4;P49019"
    uniprots[ab_targets == "NPY4R;NPY4R2"] <- "P50391;P0DQD5"
    uniprots[ab_targets == "PPAN-P2RY11;P2RY11"] <- "A0A0A6YYI3;Q96G91"
    
    # For easy matching between target and uniprot
    target_uniprot <- ab_targets %>% `names<-`(., uniprots)
    
    output$select_target <- renderUI({
      selectInput(NS(id, "select_target"), "Select Ab target", ab_targets, selected = "ADORA1")
    })
    
    output$select_uniprot <- renderUI({
      selectInput(NS(id, "select_uniprot"), "Select UniProt ID", uniprots, selected = "P30542")
    })
    
    # Sync inputs
    observeEvent(input$select_uniprot, {
      updateSelectInput(session, "select_target", "Select Ab target", ab_targets, target_uniprot[input$select_uniprot])
    })
    observeEvent(input$select_target, {
      updateSelectInput(session, "select_uniprot", "Select UniProt ID", uniprots, names(target_uniprot)[target_uniprot == input$select_target])
    })
    
    return(reactive(input$select_target))
  })
}


## Making the plots
sibcor_plot_ui <- function(id) {
  ns <- NS(id)
  plotOutput(ns("plt"), height = "700px", width = "700px")
}

sibcor_plot_server <- function(id, data_in, target_in) {
  # data_in, data containing long format data with 1D4 detection
  # target_in, reactive antibody target name (string)
  moduleServer(id, function(input, output, session) {
    
    output$plt <- renderPlot({
      # Prevent error from non-existent selection
      sel_target <- req(target_in())
      
      data_in$abval_dat %>%
        filter(gene_name == sel_target,
               gene_name_id != "CALCRL.299",
               gpcr_code != "positive IK19") %>%
        mutate(gpcr_ab = map_chr(gene_id_ab, function(x) {
          spl_name <- str_split(x, "\\.") %>% unlist()
          paste0(spl_name[1], "_", spl_name[3])
        })) %>%
        pivot_wider(id_cols = "unique_sample_name",
                    names_from = "gpcr_ab",
                    values_from = "value_rz") %>%
        column_to_rownames("unique_sample_name") %>%
        pairs.panels(method = "pearson", rug = F, pch = 21, ellipses = F, lm = T,
                     breaks = "FD", hist.col = "grey70", cex.axis = 1.5)
    }, width = 700, height = 700)
  })
}

## Antibody table
sibcor_ab_tbl_ui <- function(id) {
  ns <- NS(id)
  DTOutput(ns("ab_tbl"))
}

sibcor_ab_tbl_server <- function(id, data_in, target_in) {
  # data_in, data containing antibody information
  # target_in, the reactive name of the selected Ab target
  moduleServer(id, function(input, output, session) {
    output$ab_tbl <- renderDT({
      sel_target <- req(target_in())
      
      ab_tbl <- data_in$prot %>%
        filter(gene_name == sel_target, gene_name_id != "CALCRL.299") %>%
        select(antibody, gene_name, uniprot, prest, pass) %>%
        # Uniprot ids for multitargets
        mutate(uniprot = case_when(gene_name == "HCAR2;HCAR3" ~ "Q8TDS4;P49019",
                                   gene_name == "NPY4R;NPY4R2" ~ "P50391;P0DQD5",
                                   gene_name == "PPAN-P2RY11;P2RY11" ~ "A0A0A6YYI3;Q96G91",
                                   T ~ uniprot)) %>%
        # Remove letters from PrEST ID
        mutate(prest = str_extract(prest, "\\d+")) %>%
        # yes/no instead of true/false
        mutate(pass = case_when(pass ~ "Yes", T ~ "No")) %>%
        rename("GPCR" = "gene_name", "UniProt" = "uniprot", "Antibody" = "antibody",
               "Antigen" = "prest", "On-target" = "pass")
      
      datatable(ab_tbl, rownames = F, options = list(dom = "tr", pageLength = nrow(ab_tbl)))
    })
  })
}


### Session information ###
session_info_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("System information"),
    htmlOutput(ns("system_info")),
    h4("Added packages"),
    tableOutput(ns("other_pckg_versions"))
  )
}

session_info_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    s <- sessionInfo()
    
    output$system_info <- renderUI({
      paste0(
        s$R.version$version.string, "<br>",
        "Platform: ", s$platform, "<br>",
        "Running under: ", s$running
      ) %>% HTML()
    })
    
    # Function for making a table of packages, their versions, and the citations from the citation() function
    pckg_version_tbl <- function(pckg_names) {
      pckg_versions <- sapply(pckg_names, getNamespaceVersion)
      
      # Also get citations for each package
      pckg_citations <- sapply(pckg_names, function(x) {
        # Get citation in HTML format, edit links to open in new tabs
        citation(x) %>% format(style = "html") %>% str_replace_all("<a href", "<a target='_blank' rel='noopener noreferrer' href") %>%
          unlist() %>% paste(collapse = "")
      })
      
      tbl_out <- data.frame("Package" = pckg_names,
                            "Version" = pckg_versions,
                            "Citation" = pckg_citations)
      return(tbl_out)
    }
    
    # Other attached packages (not base)
    output$other_pckg_versions <- renderTable({
      pckg_version_tbl(sort(names(s$otherPkgs)))
    }, striped = T, sanitize.text.function = function(x) x)
    
  })
}
