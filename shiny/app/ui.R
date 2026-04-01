library(shiny)
library(bslib)
library(shinyWidgets)
library(shinycssloaders)


ps_intro <- fluidPage(
  layout_column_wrap( 
    width = 1/2,
      card(
        card_header("Introduction"),
        tags$p("This shiny app shows the main and supplementary figures for the manuscript ", 
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
        card_image("../../code_R_analysis/output_plots/fig0_shiny.svg", height = "600px"),
        card_footer("The GMGC dataset was analyzed in nucleotide and/or amino acid format by the ARG detection tools to obtain a list of putative resistance genes. We quantified the resistome (abundance, diversity, pan- and core-resistome) using the putative ARGs in host-associated and external environments.")
        )
      )
  )
  

# ARGs Tab
ps_args <- page_sidebar(
  sidebar = sidebar(
    
    pickerInput(
      inputId = "tools_unigenes",
      label = "Choose the tools you want to compare:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} tools selected"
      )
    )
    ),
  
  layout_column_wrap( 
    width = 1, fill = TRUE,
        card(
          full_screen = TRUE,  # optional but nice
          fill = TRUE,  
          card_header("(a) Number of putative ARGs and (b) gene class proportion"),
          withSpinner(plotOutput("plot_count_genes_tool", height = "600px"), type = 8, color = "#1b9e77") # Spinner + Stable height 
        )
      )
)

# Pan- and Core-resistome Tab
ps_pan_core <- page_sidebar(
  sidebar = sidebar(

    pickerInput(
      inputId = "tool_pan_core",
      label = "Choose the tools you want to show the number of genes for:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} tools selected"
      )
    ),
    
    pickerInput(
      inputId = "environment_pan_core",
      label = "Choose the environments you want to show:",
      choices = as.list(EN),
      selected = EN[c(1,9,10,13)],
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} environments selected"
      )
    ),
    
    pickerInput(
      inputId = "single_environment_pan_core",
      label = "Choose the environment you want to show the proportion of genes for:",
      choices = as.list(EN),
      selected = "human gut",
      multiple = FALSE,
      options = list(
        `live-search` = TRUE # Adds a search bar for convenience
      )
    ),
    
    pickerInput(
      inputId = "threshold_proportion",
      label = "Proportion of metagenomic samples where the gene appears in each subsampling event:",
      choices = list(">= 30%" = 0.3, ">= 40%" = 0.4, ">= 50%" = 0.5, ">= 60%" = 0.6, ">= 70%" = 0.7, ">= 80%" = 0.8, ">= 90%" = 0.9),
      selected = 0.5,
      multiple = FALSE
    ),
    
    pickerInput(
      inputId = "threshold_samples",
      label = "Minimum number of subsets the gene has to be part of the subsample core-resistome:",
      choices = list(
        ">= 200" = 200,
        ">= 250" = 250,
        ">= 300" = 300,
        ">= 350" = 350,
        ">= 400" = 400,
        ">= 450" = 450,
        "500" = 500
      ),
      selected = 450,
      multiple = FALSE
    ),
    
    radioGroupButtons(
      inputId = "threshold_pan_core_id",
      label = "Identity threshold DeepARG/RGI (amino acid):",
      choices = c(
        "Default" = "0.0",
        ">= 60%" = "60.0",
        ">= 70%" = "70.0",
        ">= 80%" = "80.0"
      ),
      selected = "0.0",
      status = "primary",
      size = "sm",
      justified = TRUE
    )
),
  
  
  layout_column_wrap( 
    width = 1,
      card(
        card_header("Number of genes"),
        withSpinner(plotOutput("pan_core", height = "600px"), type = 8, color = "#1b9e77") 
      )
  )
)


## Abundance and Diversity Tab

ps_abundance_diversity <- page_sidebar(
  sidebar = sidebar(

    pickerInput(
      inputId = "tool_abundance",
      label = "Choose the tools you want to show:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,            # Adds "Select All / Deselect All" buttons
        `selected-text-format` = "count > 2", # Shows "X items selected" if >2 are chosen
        `count-selected-text` = "{0} tools selected"
      )
    ),
    
    pickerInput(
      inputId = "environment_abundance",
      label = "Choose the environments you want to show:",
      choices = as.list(EN),
      selected = EN[1],
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} environments selected"
      )
    ),
    
    
    pickerInput(
      inputId = "abundance_genes",
      label = "Choose the genes you want to show:",
      choices = as.list(as.character(gene_classes)),
      selected = top_abundance,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,             # Adds a search bar to quickly find genes!
        `selected-text-format` = "count > 3",
        `count-selected-text` = "{0} genes selected"
      )
    ),
  
  radioGroupButtons(
    inputId = "plot_other",
    label = "Plot other gene classes together:",
    choices = c("Yes", "No"),
    selected = "Yes",
    status = "primary",
    size = "sm",
    justified = TRUE
  )
  ),
  
  navset_card_underline(
    nav_panel(
      "Relative abundance",
      fill = TRUE,
      page_fillable(
        layout_columns(
          #col_widths = c(6, 6, 12),
          #row_heights = c("1fr", "0.3fr"),
          
          card(
            fill = TRUE,
            card_header("Relaltive abundance per sample and gene class"),
            withSpinner(plotOutput("plot_abundance", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
          )
        )
      )
    ),
    
    ## Diversity menu
    nav_panel(
      "Richness",
      page_fillable(
        layout_columns(
          col_widths = c(6, 6, 12),
          row_heights = c("1fr", "0.3fr"),
          
          card(
            card_header("Richness per sample"),
            #withSpinner(plotOutput("plot_diversity", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
          )
        )
      )
    )
  )
)

## Overlaps Tab
ps_overlap <- page_sidebar(
  sidebar = sidebar(
    pickerInput(
      inputId = "tool_overlap",
      label = "Choose the tools you want to show:",
      choices = tool_choices,
      selected = c("DeepARG","fARGene", "RGI-DIAMOND"),
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,           
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} tools selected"
      )
    ),
    
    pickerInput(
      inputId = "tool_overlap_calc",
      label = "Choose the tools you want to include in the calculation:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} tools selected"
      )
    ),
    
    pickerInput(
      inputId = "overlap_genes",
      label = "Choose the genes you want to show:",
      choices = as.list(as.character(gene_classes)),
      selected = top_cso,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,             
        `selected-text-format` = "count > 3",
        `count-selected-text` = "{0} genes selected"
      )
    )
  ), 
  
  navset_card_underline(
    
    # Overview panel
    nav_panel(
      "CSC",
      layout_columns(
        #$col_widths = c(6, 6, 6, 6, 12),
        card(
          full_screen = TRUE,
          withSpinner(plotOutput("overlap", height = "450px"), type = 8, color = "#1b9e77")
        )
      )
    )
  )
)

  
# Define UI for the argCompare application
page_navbar(
  theme = "yeti",
  # theme = bs_theme(bootswatch = "yeti", primary = "#1b9e77"),
  inverse = TRUE,
  nav_panel(title = "Introduction", p(ps_intro)),
  nav_panel(title = "ARGs", p(ps_args)),
  nav_panel(title = "Abundance and diversity", p(ps_abundance_diversity)),
  nav_panel(title = "Pan- and core-resistome", p(ps_pan_core)),
  nav_panel(title = "Overlap", p(ps_overlap)),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Big Data Biology Lab", href = "https://www.big-data-biology.org")),
    nav_item(tags$a("CMR", href = "https://research.qut.edu.au/cmr/"))
  )
)
