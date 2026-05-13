responsive_css <- tags$head(
  tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0"),
  tags$style(HTML("

    html, body {
      overflow-x: hidden;
      -webkit-text-size-adjust: 100%;
    }

    .html-fill-item,
    .html-fill-container {
      min-width: 0 !important;
      max-width: 100% !important;
    }

    .card {
      width: 100%;
      min-width: 0;
      max-width: 100%;
      box-sizing: border-box;
      word-break: break-word;
      height: auto !important;
    }
    .card-body {
      width: 100% !important;
      min-width: 0 !important;
      box-sizing: border-box;
      height: auto !important;
      overflow: visible !important;
    }

    .shiny-plot-output img,
    .shiny-plot-output canvas {
      max-width: 100% !important;
    }

    .card-img, .card-img-top, .card-img-bottom {
      max-width: 100%;
      height: auto;
    }

    .bootstrap-select .dropdown-menu {
      max-width: 100%;
      min-width: 0;
      word-wrap: break-word;
      white-space: normal;
    }

    .sidebar .control-label,
    .sidebar .shiny-input-container {
      max-width: 100%;
      word-break: break-word;
    }

    .bslib-page-sidebar,
    .bslib-sidebar-layout,
    .bslib-sidebar-layout > .bslib-sidebar-main,
    .bslib-sidebar-layout > .main,
    .tab-content,
    .tab-pane,
    .bslib-page-fill {
      height: auto !important;
      max-height: none !important;
      overflow: visible !important;
    }

    .navbar-nav { flex-wrap: wrap; }

    @media (min-width: 1400px) {
      .layout-column-wrap {
        display: flex !important;
        flex-direction: row !important;
        width: 100% !important;
      }
      .layout-column-wrap > * {
        flex: 0 0 50% !important;
        width: 50% !important;
        max-width: 50% !important;
        min-width: 0 !important;
      }
      .card-img { max-height: 900px; object-fit: contain; }
    }

    @media (max-width: 1399px) {
      .layout-column-wrap {
        display: flex !important;
        flex-direction: column !important;
        width: 100% !important;
      }
      .layout-column-wrap > *,
      .layout-column-wrap > .html-fill-item {
        flex: 0 0 100% !important;
        width: 100% !important;
        max-width: 100% !important;
        min-width: 0 !important;
        margin-left: 0 !important;
        margin-right: 0 !important;
      }
      .bslib-sidebar-layout {
        flex-direction: column !important;
      }
      .bslib-sidebar-layout > .sidebar {
        width: 100% !important;
        max-width: 100% !important;
        border-right: none !important;
        border-bottom: 1px solid #dee2e6;
        padding-bottom: .75rem;
      }
      .card-img { max-height: 700px; object-fit: contain; }
    }

    @media (max-width: 991px) {
      .card-img { max-height: 600px; object-fit: contain; }
    }

    @media (max-width: 767px) {
      .container-fluid {
        padding-left: .5rem !important;
        padding-right: .5rem !important;
      }
      .navbar-nav { width: 100%; }
      .navbar-nav .nav-item { width: 100%; text-align: center; }
      .card-body   { padding: .6rem !important; }
      .card-header { padding: .5rem .75rem !important; font-size: .95rem; }
      .card-footer { padding: .5rem .75rem !important; font-size: .8rem; }
      .bootstrap-select,
      .bootstrap-select > .dropdown-toggle { width: 100% !important; }
      .sidebar em { font-size: .8rem; }
      .card-img { max-height: 420px; object-fit: contain; }
      .card-body p,
      .card-body li { font-size: .9rem; line-height: 1.5; }
    }

    @media (max-width: 479px) {
      .navbar-brand { font-size: 1rem; }
      .container-fluid {
        padding-left: .25rem !important;
        padding-right: .25rem !important;
      }
      .card-img { max-height: 320px; }
      .card-body p,
      .card-body li { font-size: .85rem; }
    }

  "))
)


ps_intro <- fluidPage(
  tags$style(HTML("
    .card > .card-header { font-size: 1.3rem !important; font-weight: 600 !important;}
    .card p { font-size: 0.95rem; line-height: 1.6; }
    .card ul { font-size: 0.95rem; line-height: 1.8; }
    .card-footer { font-size: 0.85rem; color: #666; }
  ")),
  
  layout_column_wrap( 
    width = 1/2,

    card(
      card_header(
        h3(
          "Antibiotic Resistance Gene Detection on the Global Microbial Gene Catalog v1.0 (",
          tags$a(
            "GMGC",
            href = "https://gmgc.embl.de/",
            target = "_blank"
          ),
          ") dataset"
        )
      ),
      
      card_body(
        p(
          "This dataset supports the manuscript ",
          tags$b(em("The elusive resistome: a global comparison reveals large discrepancies among detection pipelines.")),
          tags$br(),
          tags$br(),
          tags$p(
            "The preprint is available on",
            tags$a(
              "BioRxiv",
              href = "https://www.biorxiv.org/content/10.64898/2026.05.11.724158v1",
              target = "_blank"
            ),
            ". In addition to this, a full description of the datasets and access to them are available in ",
            tags$a(
              "Zenodo",
              href = "https://zenodo.org/records/19702877",
              target = "_blank"
            ),
            "."
          )
        ),
        
        tags$hr(),
        
        h4("Unigenes and metagenomes"),
        
        tags$ul(
          tags$li(
              "ARG predictions were done on 278,788,551 unigenes"
          ),
          tags$li(
            "The abundance and richness of the ARGs was estimated on 11,519 metagenomic samples"
          ),
          tags$li(
            "The pan- and core-resistomes were estimated across 13 different habitats represented by the metagenomic samples"
          )
        ),

        tags$hr(),
        
        h4("Detection Pipelines"),
        
        tags$ul(
          tags$li(
            tags$a(
              "fARGene v0.1",
              href = "https://github.com/fannyhb/fargene",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              "DeepARG v2",
              href = "https://github.com/gaarangoa/deeparg",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              "AMRFinderPlus v4.0.15",
              href = "https://github.com/ncbi/amr",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              "RGI v6.0.3 (CARD v4.0.0)",
              href = "https://github.com/arpcard/rgi",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              "ResFinder v2.4.0",
              href = "https://github.com/cadms/resfinder",
              target = "_blank"
            )
          ),
          tags$li(
            tags$a(
              "ABRicate v1.0.1",
              href = "https://github.com/tseemann/abricate",
              target = "_blank"
            )
          )
        ),
        
        tags$hr(),
        
        h4("Ontology harmonization"),
        
        tags$ul(
          tags$li(
            "Outputs from DeepARG, AMRFinderPlus, ABRicate, and ResFinder were standardized using",
            tags$a(
              "argNorm v1.0.0",
              href = "https://github.com/BigDataBiology/argNorm",
              target = "_blank"
            ),
            "."
          )
        ),
      )
      ),
        
    card(
      card_header("Pipeline"),
      card_image("../../code_R_analysis/output_plots/fig0_shiny.svg", height = "800px"),
      card_footer("We refer readers to the manuscript for details on the methods. The icons", 
      em("excess-weight-male-clothed, mouse-gray, pig-white,"), 
      "and", 
      em("dog,"),
      "in the abstract figure were adapted from",
      tags$a(
        "Servier Medical Art",
        href = "https://smart.servier.com",
        target = "_blank"
      ),
      "and are licensed under",
      tags$a(
        "CC BY 3.0 Unported",
        href = "https://creativecommons.org/licenses/by/3.0/",
        target = "_blank"
      )
    )
)
)
)

# ARGs Tab
ps_args <- page_sidebar(
  fillable = FALSE, 
  sidebar = sidebar(
    width = 300,
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
      label = "Filter gene classes (max 20):",
      choices = gene_classes,
      selected = gene_classes_default,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `max-options` = 20,
        `max-options-text` = "Limit reached (max 20 gene classes")
    ),
    
    markdown( "*Use the sidebar to filter specific pipelines or focus on particular gene classes to see how the resistome profile shifts.*")
  ),
  
  nav_panel(
    "Number of ARGs and Gene Class Proportion",
      layout_columns(
        col_widths = breakpoints(xs = 12, xxl = 6),
        
      layout_columns(
        col_widths = 12,
        
        card(
          card_header("Number of ARGs"),
          markdown(
            "This section highlights the difference in number of Antimicrobial Resistance Genes (ARGs) reported by each pipeline.
          \n
          A total of 178,107 unigenes from GMGCv1 were reported as ARG by at least one pipeline.
          The largest difference, 45-fold, was observed between ABRicate-ResFinder and DeepARG."
          ),
          withSpinner(plotOutput("plot_count_genes_tool", height = "500px"), type = 8, color = "#1b9e77"),
          downloadButton("download_count_genes", "Download Table")
        ),
        
        card(
          card_header("Jaccard Index"),
          markdown(
            "This is the Pairwise Jaccard index showing the mean overlap in ARGs detected between pipelines.
            \n
            Higher values indicate greater agreement between two pipelines in the set of genes they identify as ARGs."
          ),
          withSpinner(plotOutput("jaccard_plot", height = "500px"), type = 8, color = "#1b9e77"),
          downloadButton("download_jaccard_plot", "Download Table")
        )
      ),
      
        card(
          card_header("ARG Class Proportion"),
          markdown(
            "A heatmap illustrating the proportion of ARG classes as reported by each pipeline. Note that for 
             this heatmap, only classes representing 5% of the pipeline's total have a label. 
             \n
             This plot reveals how the proportional makeup of gene classes shifts depending on the chosen pipeline."
          ),
          withSpinner(plotOutput("plot_gene_class_proportion", height = "700px", fill = TRUE), type = 8, color = "#1b9e77"),
          downloadButton("download_gene_class_proportion", "Download Table")
        )
      )
    )
  )


# Abundance and Diversity Tab

ps_abundance <- page_sidebar(
  fillable = FALSE,
  sidebar = sidebar(
    width = 300,
    
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
      label = "Gene classes to show (max 10):",
      choices = as.list(as.character(gene_classes)),
      selected = top_abundance,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,             
        `selected-text-format` = "count > 3",
        `count-selected-text` = "{0} gene classes selected",
        `max-options` = 10,
        `max-options-text` = "Limit reached (max 10 gene classes)"
      )
    ),
    
    markdown("*Use this sidebar to filter specific pipelines, habitats, or focus on particular gene classes to observe these shifts/difference in abundance of genes estimation.*")
  ),
  

  navset_card_underline(
    nav_panel(
      "Relative Abundance per Sample and Gene Class",
      layout_columns(
        col_widths = breakpoints(xs = 12, xxl = 6),
      
          card(
              card_header("Relative abundance per Sample"),
              markdown(
                "
                We show the distribution of the relative abundance of ARGs detected by each pipeline across habitats. The middle line 
                denotes the median while each box limits represent the interquartile range and the
                whiskers extend to 1.5×IQR beyond the first and third quartiles."),
              
              withSpinner(plotOutput("plot_abundance", height = "700px"), type = 8, color = "#1b9e77"),
              downloadButton("download_abundance", "Download Table")
            ),
          
          card(
            card_header("Relative Abundance per Gene Class"),
            markdown(
              "
            We show the distribution of the relative abundance of ARG classes detected by each pipeline across habitats. The middle line 
            denotes the median while each box limits represent the interquartile range. ARG class abbreviations are found in the supplementary tables."
            ),
            withSpinner(plotOutput("plot_abundance_gene_class", height = "700px"), type = 8, color = "#1b9e77"),
            downloadButton("download_class_abundance", "Download Table")
          )
        )
    ),
    
    nav_panel(
      "Richness",
      layout_column_wrap( 
        width = 1,
        card(
          
          card_header("Richness"),
          markdown("This plot shows the **Richness** of ARGs detected by each pipeline per sample across habitats."),
          withSpinner(plotOutput("plot_richness", height = "700px"), type = 8, color = "#1b9e77"),
          downloadButton("download_richness", "Download Table")
        )
      )
    )
  )
)

# Pan- and Core-resistome Tab
ps_pan_core <- page_sidebar(
  fillable = FALSE,
  sidebar = sidebar(
    width = 300,
    
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
    
    pickerInput(
      inputId = "threshold_samples",
      label = "Minimum number of subsets the gene has to be part of the subsample core-resistome (n):",
      choices = list(
        ">= 200" = 200,
        ">= 250" = 250,
        ">= 300" = 300,
        ">= 350" = 350,
        ">= 400" = 400,
        ">= 450" = 450
      ),
      selected = 450,
      multiple = FALSE
    ),
    
    pickerInput(
      inputId = "threshold_proportion",
      label = "Proportion of metagenomic samples where the gene appears in each core-subsample (p):",
      choices = list(">= 30%" = 0.3, ">= 40%" = 0.4, ">= 50%" = 0.5, ">= 60%" = 0.6, ">= 70%" = 0.7, ">= 80%" = 0.8, ">= 90%" = 0.9),
      selected = 0.5,
      multiple = FALSE
    ),
 
    markdown("*Use the sidebar to compare specific pipelines, select different habitats, or adjust the strictness of the core-resistome threshold.*")
  ),
  
  layout_column_wrap( 
    width = 1,
    card(
      card_header("Number of genes"),
      markdown(
        "This section shows the size of the **Pan-** and **Core** resistomes across different habitats.
        
        For each pipeline and habitat, the core-resistome was estimated by randomly selecting 500 subsamples of 100 metagenomic samples. 
        For each subsample, we recorded: 
        1. the subsample core, i.e., the gene class of the centroids with a detection value ≥1 in at least *p=50%* of samples, grouped by habitat, 
        2. the total number of ARGs with a detection value ≥1 in any sample, as a measure of richness. 
        
        The pan-resistome size for each habitat and ARG class was calculated as the mean richness across the 500 subsamples. 
        The core-resistome corresponded to ARGs present in *≥n=90%* of the subsample cores.
        
        Alternative values for *p* and *n* can be selected in the menu on the left." 
      ),
      withSpinner(plotOutput("pan_core", height = "750px"), type = 8, color = "#1b9e77"),
      downloadButton("download_pan_resistome", "Download Pan-resistome Table"),
      downloadButton("download_core_resistome", "Download Core-resistome Table")
      )
    )
  )


# Overlaps Tab
ps_overlap <- page_sidebar(
  fillable = FALSE,
  sidebar = sidebar(
    width = 300,
    
    pickerInput(
      inputId = "tool_overlap",
      label = "Reference pipeline:",
      choices = tool_choices_single,
      selected = c("DeepARG", "RGI-DIAMOND", "fARGene", "AMRFinderPlus", "ResFinder"),
      multiple = TRUE,
      options = list(`actions-box` = TRUE)
    ),
    
    pickerInput(
      inputId  = "tool_overlap_comp",
      label    = "Pipelines to compare against:",
      choices  = tool_choices_single,
      selected = basic_tools,
      multiple = TRUE,                        
      options  = list(
        `actions-box`            = TRUE,
        `selected-text-format`   = "count > 2",
        `count-selected-text`    = "{0} pipelines selected"
      )
    ),
    
    pickerInput(
      inputId = "overlap_genes",
      label = "Gene classes to show (max 15):",
      choices = as.list(as.character(gene_classes)),
      selected = top_cso,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        `live-search` = TRUE,
        `selected-text-format` = "count > 3",
        `count-selected-text` = "{0} gene classes selected",
        `max-options` = 15,
        `max-options-text` = "Limit reached (max 15 gene classes)"
      )
    ),
    markdown("*Use the sidebar to select the pipelines you want to compare or filter specific gene classes to see exactly where the pipelines agree or diverge.*")
  ),
  
  nav_panel(
    "Class-Specific Coverage (CSC)",
    layout_columns(
      width  = 1,
        
        card(
          card_header("Class-Specific Coverage (CSC) by Gene Class"),
          markdown(
            "For a given ARG class, the CSC of pipeline *A* (reference pipeline) with respect to pipeline *B* is the proportion of ARGs reported by pipeline *B* that were also reported by pipeline *A*. 
            
            While analogous to recall, we intentionally use the term 'coverage' to avoid implying that any single pipeline constitutes a ground truth. 
            
            Here, we plot the distribution of CSC of a reference pipeline and gene class when compared to *n* more pipelines.
            
            A distribution closer to 100% indicates that the reference pipeline reports most of the ARGs of that class that other pipelines report."
          ),
          withSpinner(plotOutput("overlap_gene_class", height = "550px"), type = 8, color = "#1b9e77"),
          downloadButton("download_overlap_gene_class", "Download Table")
        )
      )
    )
  )

#SUPPLEMENTARY TABLES TAB
ps_tables <- page_fillable(
  navset_card_underline(
    nav_panel(
      "Sample Metadata",
      card(
        card_header("Metagenomic samples per habitat"),
        markdown("Overview of all metagenomic samples and their associated habitats."),
        reactableOutput("table_s1"),
        downloadButton("download_table_s1", "Download Table")
      )
    ),
    
    nav_panel(
      "ARG Classes",
      card(
        card_header("ARG Classes and Abbreviations"),
        markdown("ARO term IDs, gene class names, and abbreviations used throughout the app."),
        reactableOutput("table_s2"),
        downloadButton("download_table_s2", "Download Table")
      )
    ),
    
    nav_panel(
      "ARGs per Pipeline",
      card(
        card_header("ARGs Identified per Pipeline"),
        markdown("Full list of ARGs (unigenes) identified by each pipeline, with their associated gene class."),
        withSpinner(reactableOutput("table_s3"), type = 8, color = "#1b9e77"),
        downloadButton("download_table_s3", "Download Table")
      )
    )
  )
)

# Define UI for the argCompare application
page_navbar(
  theme = "yeti",
  inverse = TRUE,
  header = responsive_css,
  nav_panel(title = "Introduction", ps_intro),
  nav_panel(title = "ARGs", ps_args),
  nav_panel(title = "Abundance and Richness", ps_abundance),
  nav_panel(title = "Pan- and core-resistome", ps_pan_core),
  nav_panel(title = "Class-specific Overlap", ps_overlap),
  nav_panel(title = "Supplementary Tables", ps_tables),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Big Data Biology Lab", href = "https://www.big-data-biology.org")),
    nav_item(tags$a("CMR", href = "https://research.qut.edu.au/cmr/"))
  )
)
