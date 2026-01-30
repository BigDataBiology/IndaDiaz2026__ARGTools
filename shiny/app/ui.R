


ps_intro <- fluidPage(
  card_header(tags$b(tags$i("Introduction"))),
  full_screen = TRUE, 
  height = "100%", # Occupy full vertical space available
  

  card_body(
    fillable = TRUE, 
    padding = 0,
    layout_column_wrap(
      width = 1/2, 
      fill = TRUE,      # This forces the columns to fill the card height
      heights_equal = "row",
      
      # Sub-card for Plot A
      card(
        tags$p("This shiny app shows the main and supplementary figures form the manuscript ", 
               tags$b(tags$i("How ARG Detection Tools Shape Our View of the Resistome")), 
               "and allows to control different parameters."),
        tags$p(tags$b("302,655,267 unigenes")),
        tags$p("The unigenes are representative sequences after clustering 2.3 billion bacterial genes at 95% identity. The unigenes come from ",  tags$a(href="https://gmgc.embl.de/", "GMGC"),
               "and are accessible", tags$a(href="https://gmgc.embl.de/download.cgi", "here.")),
        tags$p(tags$b("13,174 metagenomic samples from 16 different habitats")),
        tags$p("The metagenomic samples come from GMGC and are available ",
               tags$a(href="https://gmgc.embl.de/downloads/v1.0/metadata/GMGC10.sample.meta.tsv.gz", "here"), ". In this app, we did not consider the habitats amplicon, isolate, and built-environment.",
               "The abundance of each gene in the metagenomes can be accessed", tags$a(href="https://gmgc.embl.de/downloads/v1.0/GMGC10.sample-abundance.tsv.xz", "here"), ". The summary of metagenomic samples per habitat can be found in Table S2."),
        tags$p(tags$b("ARG classes")),
        tags$p("Ontology normalization was done with", tags$a(href="https://github.com/BigDataBiology/argNorm", "argNorm"), ". Gene classes were manually curated after. The gene classes used in this project can be found in Table S3."),
        tags$p(tags$b("Tools")),
        "We used full-sized gene in all tools. For each tool, we chose a single parameter:",
        tags$ul(
          tags$li("Nucleotide sequences through:"),
          tags$ul(
            tags$li("ResFinder, and"), 
            tags$li("ABRicate with the databases: CARD, ResFinder, NCBI, ARG-ANNOT, and MEGARES 2.0.;")
          ),
          tags$li("Amino acid sequences through:"),
          tags$ul(
            tags$li("DeepARG"),
            tags$li("fARGene"), 
            tags$li("RGI with DIAMOND aligner"), 
            tags$li("AMRFinderPlus.")
          )
        )
      ),
      
      card(
        card_header("Pipeline"),
        card_image("../../code_R_analysis/output_plots/fig1.svg", height = "100%"),
        card_footer("The GMGC dataset was analyzed in nucleotide and/or amino acid format by the ARG detection tools to obtain a list of putative resistance genes. We quantified the resistome (abundance, diversity, pan- and core-resistome) using the putative ARGs in host-associated and external environments.")
        )
      )#,
    #card_footer("Note:")
    )
  )


ps_args <- page_sidebar(
  #title = "Genes detected as ARGS",
  sidebar = sidebar(
    helpText(
      ""
    ),
    selectInput(
      "tools_unigenes",
      "Choose the tools you want to compare:",
      tool_choices,
      as.vector(tools_levels),
      multiple = TRUE
    )
  ),
  card(
    full_screen = TRUE, 
    height = "100%", # Occupy full vertical space available
    card_header("Inter-tool differences in the number and classes of ARGs detected."),
    
    # card_body with fill = TRUE allows the contents to expand dynamically
    card_body(
      fillable = TRUE, 
      padding = 0,
      layout_column_wrap(
        width = 1/2, 
        fill = TRUE,      # This forces the columns to fill the card height
        heights_equal = "row",
        
        # Sub-card for Plot A
        card(
          card_header("A: Number of putative ARGs"),
          plotOutput("plot_count_genes_tool", height = "100%") # Height 100% is key
        ),
        
        # Sub-card for Plot B
        card(
          card_header("B: Proportion of ARGs per class"),
          plotOutput("plot_alluvial_classes", height = "100%") # Height 100% is key
        )
      )
    ),
    
    card_footer("Note: In B, genes forming at least 99% and with 0.5% proportion are shown.")
  )
  
)


ps_pan_core <- page_sidebar(
  #title = "Pan and Core resistome",
  sidebar = sidebar(
    helpText(
      ""
    ),
    selectInput(
      "tool_pan_core",
      "Choose the tools you want to compare:",
      tool_choices,
      as.vector(tools_levels),
      multiple = TRUE
    ),
    selectInput(
      "environment_pan_core",
      "Choose the tools you want to compare:",
      as.list(EN2),
      selected = EN2[c(1,9,10,13)],
      multiple = TRUE
    ),
    radioButtons(
      "threshold_proportion",
      "Proportion of metagenomic samples where the gene appears in each subsampling event",
      choices = list(">= 30%" = 0.3, ">= 40%" = 0.4, ">= 50%" = 0.5, ">= 60%" = 0.6, ">= 70%" = 0.7, ">= 80%" = 0.8, ">= 90%" = 0.9),
      selected = 0.5
    ),
    sliderInput(
      "threshold_samples",
      "Minimum number of subsets the gene has to be part of the subsample core-resistome",
      min = 200,
      max = 499,
      value = 450
    )
  ),
  
  card(
    full_screen = TRUE, 
    height = "100%", # Occupy full vertical space available
    card_header("Pan- and core-resitomes."),
    card_body(
      fillable = TRUE, 
      padding = 0,
      layout_column_wrap(
        width = 1, 
        fill = TRUE,
        heights_equal = "row",
        card(
          card_header("Again pan- core-resistomes"),
          plotOutput("plot_pan_core_resistome", height = "100%") # Height 100% is key
        )
      )
    ),
    
    card_footer("Note:")
  )
)



# Define UI for the argCompare application
page_navbar(
  theme = "yeti",
  #title = "How ARG Detection Tools Shape Our View of the Resistome",
  #bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "Introduction", p(ps_intro)),
  nav_panel(title = "ARGs", p(ps_args)),
  nav_panel(title = "Pan-/core-resistome", p(ps_pan_core)),
  # nav_panel(title = "Core-/pan-resistome", p(psb3)),
  # nav_panel(title = "Overlap", p(psb4)),
  # nav_panel(title = "Table S1", p(tab1)),
  # nav_panel(title = "Table S2", p(tab2)),
  # nav_panel(title = "Table S3", p(tab3)),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Big Data Biology Lab", href = "https://www.big-data-biology.org")),
    nav_item(tags$a("CMR", href = "https://research.qut.edu.au/cmr/"))
  )
)



# Define server logic 
# navset_card_underline(
#   title = "Inter-tool differences in the number and classes of ARGs detected.",
#   nav_panel("A", "The number of putative ARGs detected by each tool.", 
#             plotOutput("plot_count_genes_tool")),
#   nav_panel("B", "Proportion of ARGs per class identified by each tool. Genes forming at least 99% and with 0.5% proportion are shown.", 
#             plotOutput("plot_alluvial_classes"))
# )
