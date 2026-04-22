library(shiny)
library(bslib)
library(shinyWidgets)
library(shinycssloaders)


ps_intro <- fluidPage(
  layout_column_wrap( 
    width = 1/2,
    card(
      card_header("Introduction"),
      tags$p("This app shows the main and supplementary figures for the manuscript and allows to create different variants",
             tags$b(tags$i("How ARG Detection Pipelines Shape Our View of the Resistome")), 
             "and allows to control different parameters."),
      tags$p(tags$b("302,655,267 unigenes")),
      tags$p("The unigenes are representative sequences after clustering 2.3 billion bacterial genes at 95% identity. The unigenes come from ",  tags$a(href="https://gmgc.embl.de/", "GMGC"),
             "and are accessible", tags$a(href="https://gmgc.embl.de/download.cgi", "here.")),
      tags$p(tags$b("13,174 metagenomic samples from 16 different habitats")),
      tags$p("The metagenomic samples come from GMGC and are available ",
             tags$a(href="https://gmgc.embl.de/downloads/v1.0/metadata/GMGC10.sample.meta.tsv.gz", "here"), ". In this app, we did not consider the habitats amplicon, isolate, and built-environment.",
             "The abundance of each gene in the metagenomes can be accessed", tags$a(href="https://gmgc.embl.de/downloads/v1.0/GMGC10.sample-abundance.tsv.xz", "here"), ". The summary of metagenomic samples per habitat can be found in Table S2 in the."),
      tags$p(tags$b("ARG classes")),
      tags$p("Ontology normalization was done with", tags$a(href="https://github.com/BigDataBiology/argNorm", "argNorm"), ". Gene classes were manually curated after. The gene classes used in this project can be found in Table S3 in the paper."),
      tags$p(tags$b("Pipelines")),
      "We used full-sized gene in all pipelines. For each pipeline, we chose a single parameter:",
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
      card_image("../../code_R_analysis/output_plots/fig0_shiny.svg", height = "800px"),
      card_footer("The GMGC dataset was analyzed in nucleotide and/or amino acid format by the ARG detection pipelines to obtain a list of putative resistance genes. We quantified the resistome (abundance, diversity, pan- and core-resistome) using the putative ARGs in host-associated and external habitats")
    )
  )
)


# ARGs Tab
ps_args <- page_sidebar(
  sidebar = sidebar(
    width = 400,
    pickerInput(
      inputId = "tools_unigenes",
      label = "Pipelines to compare:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} pipelines selected"
      )
    ),
    
    pickerInput(
      inputId = "gene_classes_filter",
      label = "Filter gene classes:",
      choices = gene_classes,
      selected = gene_classes,
      multiple = TRUE,
      options = list(`actions-box` = TRUE)
    ),
    
    markdown(
      "*Use the sidebar to filter specific pipelines or focus on particular gene classes to see how the resistome profile shifts.*"
    ),
  ),
  
  nav_panel(
    "Number of ARGs and Gene Class Proportion",
    
    page_fillable(
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_header("Number of ARGs"),
          markdown(
            "This section highlights how different computational pipelines influence the detection of Antimicrobial Resistance Genes (ARGs) from the unigene dataset. 
          \n
          A total of 178,107 unigenes from GMGCv1 were reported as antibiotic resistance genes (ARG) by at least one pipeline. 
          The total ARG count varies across pipelines; for example, in the case of ABRicate-ResFinder and DeepARG, there is about a 45-fold difference in the count."
          ),
          withSpinner(plotOutput("plot_count_genes_tool", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
        ),
        card(
          card_header("ARG Class Proportion"),
          markdown(
            "A heatmap illustrating the proportion of ARG classes as reported by each pipeline. Also note that for 
             this heatmap, only classes representing 5% of the total within at least one pipeline is shown. Despite the massive differences in absolute counts seen on the left, this plot reveals how the proportional 
             makeup of gene classes shifts depending on the pipeline chosen.
             \n
             Note: we merged MFS efflux pumps with other efflux pumps."
          ),
          withSpinner(plotOutput("plot_gene_class_proportion", height = "900px", fill = TRUE), type = 8, color = "#1b9e77")
        )
      )
    ),
  )
)


## Abundance and Diversity Tab

ps_abundance <- page_sidebar(
  sidebar = sidebar(
    width = 400,
    
    pickerInput(
      inputId = "tool_abundance",
      label = "Pipelines to show:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2", 
        `count-selected-text` = "{0} pipelines selected"
      )
    ),
    
    pickerInput(
      inputId = "environment_abundance",
      label = "Habitats to show:",
      choices = as.list(EN),
      selected = EN[1],
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} habitats selected"
      )
    ),
    
    pickerInput(
      inputId = "abundance_genes",
      label = "Gene classes to show:",
      choices = as.list(as.character(gene_classes)),
      selected = top_abundance,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,             # Adds a search bar to quickly find genes!
        `selected-text-format` = "count > 3",
        `count-selected-text` = "{0} gene classes selected"
      )
    ),
    
    markdown(
      "*Use this sidebar to filter specific pipelines, habitats, or focus on particular gene classes to observe these shifts/difference in abundance of genes estimation.*"
    ),
  ),
  
  nav_panel(
    "Relative Abundance per Sample and Gene Class",
    page_fillable(
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_header("Relative abundance per Sample"),
          markdown(
            "This section shows the **Relative Abundance** of Antimicrobial Resistance Genes (ARGs) across different habitats and pipelines.
            \n
            This displays the total relative abundance of ARGs detected by each pipeline per sample across various habitats. For interpreting the boxplot, 
            it would be helpful to note that the center line denotes the median while each box limits is the interquartile range (IQR) and the
            whiskers extent to 1.5× IQR beyond the first and third quartiles."),
          
          withSpinner(plotOutput("plot_abundance", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
        ),
        
        card(
          card_header("Relative Abundance per Gene Class"),
          markdown(
            "This plots illustrates the relative abundance by gene classes based on their specific resistance mechanisms in different habitats.
            For the gene classes, the groups 'class A' and 'tet RPG' represent Class A β-lactamases and tetracycline ribosomal protection genes, respectively.
            \n
            *From the source file, we merged MFS efflux pumps with all other efflux pumps.*"
          ),
          withSpinner(plotOutput("plot_abundance_gene_class", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
        )
      )
    ),
    
    )
  )


# Pan- and Core-resistome Tab
ps_pan_core <- page_sidebar(
  sidebar = sidebar(
    width = 400,
    
    pickerInput(
      inputId = "tool_pan_core",
      label = "Pipelines to show:",
      choices = tool_choices,
      selected = basic_tools,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} pipelines selected"
      )
    ),
    
    pickerInput(
      inputId = "environment_pan_core",
      label = "Habitats to show:",
      choices = as.list(EN),
      selected = EN[c(1,9,10,13)],
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} habitats selected"
      )
    ),
    markdown(
      "*Use the sidebar to compare specific pipelines, select different habitats, or adjust the core-resistome threshold.*"
    ),
  ),
  
  layout_column_wrap( 
    width = 1,
    card(
      card_header("Number of genes"),
      markdown(
        "This section shows the **Pan-** and **Core** resistomes across different host habitats for each pipeline.
        * **Pan-resistome (Left, Panel a):** This represents the *total pool* of unique ARGs found across all samples within a given habitat by the pipelines. It is also shown that these numbers vary depending on the pipelines use.
        It also represents the diversity of resistance genes in the chosen habitat.
        \n
        * **Core-resistome (Right, Panel b):** This represents the ARGs that are persistently found across almost *all* samples within a habitat. Pipeline-specific detection influences which genes are considered ubiquitous."
      ),
      withSpinner(plotOutput("pan_core", height = "800px"), type = 8, color = "#1b9e77")
      )
    ),
  )


## Overlaps Tab
ps_overlap <- page_sidebar(
  sidebar = sidebar(
    width = 400,
    pickerInput(
      inputId = "tool_overlap",
      label = "Pipelines to show:",
      choices = tool_choices,
      selected = c("DeepARG","fARGene", "RGI-DIAMOND"),
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,           
        `selected-text-format` = "count > 2",
        `count-selected-text` = "{0} pipelines selected"
      )
    ),
    
    pickerInput(
      inputId = "overlap_genes",
      label = "Gene classes to show:",
      choices = as.list(as.character(gene_classes)),
      selected = top_cso,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `selected-text-format` = "count > 3",
        `count-selected-text` = "{0} gene classes selected"
      )
    ),
    markdown(
      "*Use the sidebar to select the pipelines you want to compare or filter specific gene classes to see exactly where the pipelines agree or diverge.*"
    )
  ),
  
  nav_panel(
    "Class-Specific Coverage (CSC)",
    page_fillable(
      layout_columns(
        col_widths = c(6, 6),
        
        card(
          card_header("Class-Specific Coverage by Pipeline"),
          markdown(
            "This tab examines **Class-Specific Coverage (CSC)**, showing the degree of overlap or agreement between different detection pipelines. In simple terms, this is it: *For a given ARG class g, the CSC of pipeline A with respect to 
              pipeline B represents the proportion of ARGs identified by pipeline B that were also captured by pipeline A.*
              
              \n
              This is the overall percentage of ARGs detected by a reference pipeline that were also successfully identified by the compared pipeline. A higher percentage indicates strong agreement between the pipelines."
          ),
          withSpinner(plotOutput("overlap", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
        ),
        
        card(
          card_header("Class-Specific Coverage (CSC) by Gene Class"),
          markdown(
            "This plot breaks down this overlap by gene classes based on specific resistance mechanisms. This reveals whether two pipelines might agree perfectly on certain gene classes (like tetracycline resistance) but completely miss each other on others."
          ),
          withSpinner(plotOutput("overlap_gene_class", height = "600px", fill = TRUE), type = 8, color = "#1b9e77")
        )
      ),
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
  nav_panel(title = "Abundance", p(ps_abundance)),
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
