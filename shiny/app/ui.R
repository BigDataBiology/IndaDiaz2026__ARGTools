


ps_intro <- fluidPage(
  layout_column_wrap( 
    width = 1/2,
      card(
        card_header("Introduction"),
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
        card_image("../../code_R_analysis/output_plots/fig1_shiny.svg", height = "100%"),
        card_footer("The GMGC dataset was analyzed in nucleotide and/or amino acid format by the ARG detection tools to obtain a list of putative resistance genes. We quantified the resistome (abundance, diversity, pan- and core-resistome) using the putative ARGs in host-associated and external environments.")
        )
      )
  )
  


ps_args <- page_sidebar(
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
    ),
    radioButtons(
      "threshold_unigenes_id",
      "Identity threshold DeepARG/RGI (amino acid)",
      choices = list("Default" = 0.0, ">= 60%" = 60.0, ">= 70%" = 70.0, ">= 80%" = 80.0),
      selected = 0.0
    )
  ),
  layout_column_wrap( 
    width = 1/2,
        card(
          card_header("A: Number of putative ARGs"),
          plotOutput("plot_count_genes_tool", height = "100%") 
        ),
        
        card(
          card_header("B: Proportion of ARGs per class"),
          plotOutput("plot_alluvial_classes", height = "100%"), 
          card_footer("Genes forming at least 99% and with 0.5% proportion are shown.")
        )
      )
)


ps_pan_core <- page_sidebar(
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
    ),
    radioButtons(
      "threshold_pan_core_id",
      "Identity threshold DeepARG/RGI (amino acid)",
      choices = list("Default" = 0.0, ">= 60%" = 60.0, ">= 70%" = 70.0, ">= 80%" = 80.0),
      selected = 0.0
    )
  ),
  
  layout_column_wrap( 
    width = 1/2,
      card(
        card_header("Number of genes"),
        plotOutput("plot_pan_core_resistome", width = "100%") 
      ),
      card(
        card_header("Proportion of genes"),
        
      )
    )
)

ps_abundance_diversity <- page_sidebar(
  sidebar = sidebar(
    helpText(
      ""
    ),
    selectInput(
      "tool_abundance",
      "Choose the tools you want to compare:",
      tool_choices,
      as.vector(tools_levels),
      multiple = TRUE
    ),
    selectInput(
      "environment_abundance",
      "Choose the tools you want to show:",
      as.list(EN2),
      selected = EN2[c(1)],
      multiple = TRUE
    ),
    radioButtons(
      "threshold_abundance_id",
      "Identity threshold DeepARG/RGI (amino acid)",
      choices = list("Default" = 0.0, ">= 60%" = 60.0, ">= 70%" = 70.0, ">= 80%" = 80.0),
      selected = 0.0
    ),
    selectInput(
      "abundance_genes",
      "Choose the genes you want to show:",
      as.list(as.character(gene_classes)),
      top20,
      multiple = TRUE
    ),
    radioButtons(
      "plot_other",
      "Plot other gene classes together:",
      as.list(c("Yes","No")),
      "Yes"
    )
  ),
  

    page_fillable(
      #title = "Abundance and diversity",
      layout_columns(
        col_widths = c(6, 6, 6, 6, 12),
        row_heights = c("1fr", "1fr", "0.3fr"),
        card(
          card_header("Abundance per sample"),
          #"This is outcome",
          plotOutput("plot_abundance", height = "100%") 
        ),
        card(
          card_header("Diversity per sample"),
          #"This is outcome 2",
          plotOutput("plot_diversity", height = "100%") 
        ),
        card(
          card_header("Median abundance per class"),
          #"This is outcome",
          plotOutput("plot_abundance_class", height = "100%") 
        ),
        card(
          card_header("Median diversity per class"),
          plotOutput("plot_diversity_class", height = "100%") 
        ),
        card(
        #  card_header("Median diversity per class"),
          plotOutput("plot_diversity_class_legend", height = "100%") 
        )
      )
    )
)



ps_overlap <- page_sidebar(
  sidebar = sidebar(
    helpText(
      ""
    ),
    selectInput(
      "tool_overlap",
      "Choose the tools you want to show:",
      tool_choices,
      as.vector(tools_levels[c(1,2,4,5)]),
      multiple = TRUE
    ),
    selectInput(
      "tool_overlap_calc",
      "Choose the tools you want to include in the calculation:",
      tool_choices,
      as.vector(tools_levels),
      multiple = TRUE
    ),
    selectInput(
      "overlap_genes",
      "Choose the genes you want to show:",
      as.list(as.character(gene_classes)),
      top20,
      multiple = TRUE
    ),
    radioButtons(
      "threshold_overlap_id",
      "Identity threshold DeepARG/RGI (amino acid)",
      choices = list("Default" = 0.0, ">= 60%" = 60.0, ">= 70%" = 70.0, ">= 80%" = 80.0),
      selected = 0.0
    )
  ),
  
    page_fillable(
      layout_columns(
        col_widths = c(10, 2, 10, 2, 10), 
        row_heights = c("1fr", "1fr", "0.3fr"),
        card(
          card_header("CSTC"),
          plotOutput("overlap_cstc", height = "100%") 
        ),
        card(
          card_header("Medians per class"),
          plotOutput("overlap_cstc_summary", height = "100%") 
        ),
        card(
          card_header("CSNO"),
          plotOutput("overlap_csno", height = "100%") 
        ),
        card(
          card_header("Medians per class"),
          plotOutput("overlap_csno_summary", height = "100%") 
        ),
        card(
          card_header(""),
          plotOutput("overlap_legend", height = "100%") 
        )
      )
    )#,
    
    #card_footer("Note:")
  
)


# Define UI for the argCompare application
page_navbar(
  theme = "yeti",
  #title = "How ARG Detection Tools Shape Our View of the Resistome",
  #bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "Introduction", p(ps_intro)),
  nav_panel(title = "ARGs", p(ps_args)),
  nav_panel(title = "Abundance and diversity", p(ps_abundance_diversity)),
  nav_panel(title = "Pan- and core-resistome", p(ps_pan_core)),
  nav_panel(title = "Overlap", p(ps_overlap)),
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



