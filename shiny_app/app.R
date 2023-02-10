### --- GPCR Antibody validation Shiny app --- ###
## Leo Dahl


### Packages ###
library(shiny)
library(tidyverse)
library(plotly)
library(ggbeeswarm)
library(patchwork)
library(data.table)
library(DT)
library(ComplexUpset)
library(psych)
library(viridis)
library(RColorBrewer)
library(bslib)
library(sass)

### Modules ###
source("modules.R")

### Load Raw data, Zscores and robust Zscores ###
dat_list <- readRDS("gpcr_mfi_data.rds")

### UI ###
ui <- fluidPage(
  
  # Different tabs and sidepanels for different analyses
  navbarPage(
    theme = bs_theme(version = 3, bootswatch = "flatly") %>%
      # Add extra formatting for GPCR info table and
      bs_add_rules(sass_file("www/shiny_app.scss")),
    collapsible = T,
    title = NULL,
    
    tabPanel(
      strong("GPCR Antibody Validation Study"),
      
      HTML("Companion app to<h2><i>Multiplexed selectivity screening of anti-GPCR antibodies</i></h2>"),
      HTML("Leo Dahl, Ilana B. Kotliar, Annika Bendes, Tea Dodig-Crnković, Samuel Fromm, Arne Elofsson, Mathias Uhlén, Thomas P. Sakmar and Jochen M. Schwenk"),
      h3("About"),
      HTML("<p><strong>G protein-coupled receptors (GPCRs) </strong> control critical cellular signaling pathways. Therapeutic agents, such 
           as antibodies (Abs), are being developed to modulate GPCR function. However, validating the selectivity 
           of anti-GPCR Abs is challenging due to sequence similarities of individual receptors within GPCR subfamilies. 
           To address this, we developed a multiplexed immunoassay to test >400 anti-GPCR Abs from the Human Protein Atlas 
           targeting a customized library of 215 expressed and solubilized GPCRs representing all GPCR subfamilies. We found 
           that ~61% of Abs were selective for their intended target, ~11% to bind off-target, and ~28% not to bind any GPCR. 
           Antigens of on-target Abs were, on average, significantly longer, more disordered, and less likely to be buried in 
           the interior of the GPCR protein than the other Abs. These results provide important insights into the immunogenicity 
           of GPCR epitopes and form a basis for the design of therapeutic Abs and the detection of pathological auto-antibodies.</p>"),
      HTML("<p><strong>This app</strong> provides complementary visualizations for the Ab validation, such as data distributions for each individual GPCR,
           plots summarizing the on-target and off-target statuses of Abs, and comparisons between paired Abs binding the same GPCR.</p>"),
      HTML("<p><strong>Data</strong> used for the study can be found at the 
           <a href='https://figshare.scilifelab.se' target='_blank' rel='noopener noreferrer'><strong> SciLifeLab Data Repository on FigShare </strong></a>
           with the ID [Add ID!]<br>
           and <strong> code </strong> used is found at the 
           <a href='https://github.com/Schwenk-Lab/' target='_blank' rel='noopener noreferrer'><strong> Schwenk Lab </strong></a> and
           <a href='https://github.com/ElofssonLab/' target='_blank' rel='noopener noreferrer'><strong> Elofsson Lab </strong></a> GitHub accounts. </p>"),
      hr(),
      HTML("App developed by Leo Dahl"), br(),
      HTML("App version 1.0.0"), br(),
      HTML("Using the bslib Flatly theme")
    ),
  
    tabPanel(
      "Antibody validation",
      sidebarLayout(
        sidebarPanel(
          
          # Select GPCR for plotting
          gpcr_select_ui("abval_gpcr_select"),
          # Description
          HTML("<hr class='sidebarhr'>
      <details><summary><strong>Click for description</strong></summary>
      <p>
      This page provides information about each included GPCR, with links to resources like the Human Protein Atlas (HPA), Ensembl, and Uniprot.<br><br>
      <strong>Expression</strong> of each GPCR was evaluated by capturing the FLAG tag and detecting the 1D4 tag on the GPCRs. The log2 of the signal (x-axis in the plot) was 
      compared to samples without GPCR (mock) using ANOVA followed by Dunnett's multiple comparison test, classifying the expression as uncertain if the 
      p-value was higher than 0.05.<br><br>
      <strong>Validation of antibodies</strong> against each GPCR was performed by capturing the GPCR using the listed antibody (identifier in the form 
      GPCRname_Abname) and detecting the 1D4 tag on the GPCR. Related GPCRs were assayed using the same antibodies to map out
      off-target binding of related targets. The mean signal of the samples containing the target GPCR was compared to a population density-based cutoff,
      classifying each antibody as validated if
      <ul>
        <li>the mean signal was higher than the cutoff, and</li>
        <li>there was no off-target binding.</li>
      </ul>
      The cutoff plots can be customised slightly by changing point size. The table below contains information about each antibody, such as identifier
      (column name: 'Antibody'), its target GPCR ('GPCR'), the Research Resource Identifier ('RRID') of the Ab, whether it recognizes only its target ('On-target'), whether there
      is any off-target binding ('Crossreacts'), and if so, the names of the GPCRs in the crossreacting samples and their interactors ('Other target').
      </p>
      </details>")
        ),
        
        mainPanel(
          gpcr_info_ui("abval_gpcr_tbl"),
          hr(),
          expr_plot_ui("abval_expr_plot"),
          hr(),
          abval_beeswarm_ui("abval_beeswarm_main"),
          ab_info_ui("abval_ab_info")
        )
      )
    ),
    
    tabPanel(
      "AbVal UpSet plot",
      sidebarLayout(
        sidebarPanel(
          HTML("<p>
       <strong>Description</strong><br>
       <strong>The UpSet plot</strong> shows the number of anti-GPCR antibodies that fulfill different intersections of the
       categories presented in the matrix below the barplot. The categories include whether the target GPCR showed sufficient
       evidence of expression ('GPCR expression', see the 'Antibody Validation' tab), if the antibody detects its intended target
       ('On-target detection'), and if the antibody binds any GPCRs other than its target ('Crossreactivity'). The bars are colored
       green if, based on the categories, there is support for only on-target binding of the antibody, and black if it is uncertain
       if the antibody binds only its intended target. <br><br>
       <strong>The table</strong> below shows the same information as the table in 'Antibody validation', but has all the assayed antibodies included.
       The search bar can be used for filtering specific antibodies or GPCRs. 
       </p>
       ")
        ),
        
        mainPanel(
          abval_upset_ui("abval_upset"),
          ab_info_ui("upset_ab_info")
        )
      )
    ),
    
    tabPanel(
      "Median signal heatmaps",
      sidebarLayout(
        sidebarPanel(
          
          # Select subfamily/family
          fam_subfam_select_ui("med_hm_select"),
          
          # Heatmap layout
          hm_custom_ui("med_hm_layout"),
          
          # Description
          HTML("<hr class='sidebarhr'>
      <details><summary><strong>Click for description</strong></summary>
      <p>
      <strong>The interactive heatmap</strong> shows the median signal (as Z-score) of each available combination of GPCR and capture antibody,
      detected using the 1D4 tag. A GPCR and the antibodies targeting the GPCR roughly follow the diagonal (offset by GPCRs with multiple
      antibodies targeting them). Bright spots deviating from this diagonal show off-target binding. Hover with the mouse pointer on cells to
      display the GPCR, antibody, and associated signal. The heatmap contents can be changed by selecting different subfamilies or GPCR ligand
      families. The layout can be customised with the controls above.
      </p>
      </details>")
        ),
        
        mainPanel(
          med_hm_ui("med_hm_plot")
        )
      )
    ),
    
    tabPanel(
      "Paired Ab correlation",
      sidebarLayout(
        sidebarPanel(
          # Selecting Ab target for sibling correlation
          target_select_ui("sibcor_target_select"),
          
          # Description
          HTML("<hr class='sidebarhr'>
      <details><summary><strong>Click for description</strong></summary>
      <p>
      <strong>The plot</strong> shows correlation between pairs of antibodies that share the same target. The upper right triangle shows
      the Pearson correlation between the antibodies and the lower left triangle shows scatter plots of the robust Z-scores of the antibody pairs. Along
      the diagonal there are robust Z-score histograms for each antibody.<br><br>
      <strong>The table</strong> contains information about the plotted antibodies, with the columns containing the identifier of the antibody
      (Column name: 'Antibody'), the name of the target GPCR ('GPCR') and its UniProt ID ('UniProt'), an identifier of the antigen ('Antigen ID'),
      and whether the antibody recognizes only its intended target ('On-target').
      </p>
      </details>")
        ),
        
        mainPanel(
          sibcor_plot_ui("sibcor_plot"),
          sibcor_ab_tbl_ui("sibcor_tbl")
        )
      )
    ),
    
    tabPanel(
      "Session information",
      session_info_ui("session_info")
    )
  )
)


### Server ###
server <- function(input, output, session) {
  ## Antibody validation ##
  # GPCR selection
  abval_gpcr <- gpcr_select_server("abval_gpcr_select", dat_list)
  
  # Ab-val beeswarm
  # bsw_param <- abval_beeswarm_side_server("abval_beeswarm_side")
  abval_beeswarm_server("abval_beeswarm_main", dat_list, abval_gpcr)
  
  # Antibody information table
  ab_info_server("abval_ab_info", dat_list, abval_gpcr)
  
  # Expression plot for target GPCR of selected antibody
  expr_plot_server("abval_expr_plot", dat_list, abval_gpcr)
  
  # GPCR information table
  gpcr_info_server("abval_gpcr_tbl", dat_list, abval_gpcr)
  
  ## AbVal UpSet plot ##
  
  # Make UpSet plot
  abval_upset_server("abval_upset", dat_list, c("Crossreactivity", "On-target detection", "GPCR expression"))
  
  # Make table for the validation of the antibodies
  ab_info_server("upset_ab_info", dat_list, "all")
  
  ## Median signal heatmaps ##
  
  # Reactive selection of family or subfamily
  med_hm_fam <- fam_subfam_select_server("med_hm_select", dat_list)
  
  # Parameters for layout
  med_hm_layout <- hm_custom_server("med_hm_layout")
  
  # Make the heatmap
  med_hm_server("med_hm_plot", med_hm_fam[[1]], med_hm_fam[[2]], dat_list, med_hm_layout)
  
  ## Sibling Ab correlation ##
  
  # Selecting Ab target to display
  sibcor_target <- target_select_server("sibcor_target_select", dat_list)
  
  # Render the correlation plots
  sibcor_plot_server("sibcor_plot", dat_list, sibcor_target)
  
  # Make a table of the sibling correlation Abs
  sibcor_ab_tbl_server("sibcor_tbl", dat_list, sibcor_target)
  
  ### Session information ###
  session_info_server("session_info")
}


shinyApp(ui = ui, server = server)

