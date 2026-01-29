


ps_intro <- fluidPage(
  title = "How ARG Detection Tools Shape Our View of the Resistome",
  card(
    #card_header(""),
    card_image("../../code_R_analysis/output_plots/fig1.svg", height = "600px"),
    card_footer("The GMGC dataset was analyzed in nucleotide and/or amino acid format by the ARG detection tools to obtain a list of putative resistance genes. We quantified the resistome (abundance, diversity, pan- and core-resistome) using the putative ARGs in host-associated and external environments.")
  ),
  card(
    #card_header("Summary"),
    tags$p("This shiny app shows the main and supplementary figures form the manuscript ", 
           tags$b(tags$i("How ARG Detection Tools Shape Our View of the Resistome")), 
           "and allows to control different parameters."),
    tags$p(tags$b("Genes")),
    tags$p("The genes used correspond to the unigenes from",  tags$a(href="https://gmgc.embl.de/", "GMGC"),
           "and are accessible", tags$a(href="https://gmgc.embl.de/download.cgi", "here.")),
    tags$p(tags$b("Metagenomes")),
    tags$p("The table with the metagenomes and their accession numbers can be found in Table S2 and ",
           tags$a(href="https://gmgc.embl.de/downloads/v1.0/metadata/GMGC10.sample.meta.tsv.gz", "here"), ". We did not consider the environments amplicon, isolate, and built-environment for this study.",
           "The abundance of each gene in the metagenomes can be accessed", tags$a(href="https://gmgc.embl.de/downloads/v1.0/GMGC10.sample-abundance.tsv.xz", "here. The summary of metagenomic samples per habitat can be found in Table S2.")),
    tags$p(tags$b("Ontology")),
    "We provide the table that links each ARO to the gene classes used in this project Table S3.",
    tags$p(tags$b("Tools")),
    "To avoid redundance, we show the results of running the genes in:",
    tags$ul(
      tags$li("nucleotide format through:"),
      tags$ul(
        tags$li("DeepARG"),
        tags$li("fARGene"), 
        tags$li("RGI with DIAMOND aligner"), 
        tags$li("ResFinder, and"), 
        tags$li("ABRicate with the databases: CARD, ResFinder, NCBI, ARG-ANNOT, and MEGARES 2.0.;")
      ),
      tags$li("amino acid format through:"),
      tags$ul(
        tags$li("AMRFinderPlus.")
      ))
  ),
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
  # navset_card_underline(
  #   title = "Inter-tool differences in the number and classes of ARGs detected.",
  #   nav_panel("A", "The number of putative ARGs detected by each tool.", 
  #             plotOutput("plot_count_genes_tool")),
  #   nav_panel("B", "Proportion of ARGs per class identified by each tool. Genes forming at least 99% and with 0.5% proportion are shown.", 
  #             plotOutput("plot_alluvial_classes"))
  # )
  
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




# Define UI for the argCompare application
page_navbar(
  theme = "yeti",
  #title = "How ARG Detection Tools Shape Our View of the Resistome",
  #bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "Introduction", p(ps_intro)),
  nav_panel(title = "ARGs", p(ps_args)),
  # nav_panel(title = "Abundance/Diversity", p(psb2)),
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
