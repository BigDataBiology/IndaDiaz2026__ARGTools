library(shiny)
library(shinydashboard)
library(tidyverse)
library(shinyWidgets)
library(plotly)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(eulerr)
library(reactable)
library(Cairo)
library(ggalluvial)
library(cowplot)
library(scales)


options(dplyr.summarise.inform = FALSE)

source("code_R_analysis/helper.R")

# Defining constants
general_size <- 10

# Label formatting function
lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}

# Color palettes
pal_7 <- brewer.pal(8, "BrBG")
pal_7 <- pal_7[c(4,5,3,6,2,7,1)]
pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]

pal_10_complete <- brewer.pal(10, "BrBG")
pal_10_complete <- pal_10_complete[c(-1,-10)]

# Pattern settings for textures
pattern_density <- 0.001 
pattern_spacing <- 0.025
pattern_fill <- "white"
pattern_size <- 0.12

# All Habitats
EN <- c("human gut", "human oral", "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil", "amplicon", "isolate", "built-environment")

#Source for each habitat
SO <- c(rep("humans", 5), rep("mammals", 4),  
        "wastewater", "marine", "freshwater", 
        "soil", rep("other", 2), "built-environment")

names(SO) <- EN

# Tool definitions
tools_levels <- c("DeepARG", "fARGene",
                  "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                  "RGI-DIAMOND", "ABRicate-CARD",
                  "AMRFinderPlus", "ABRicate-NCBI",
                  "ResFinder", "ABRicate-ResFinder")

tools_labels <- c("DeepARG", "fARGene",
                  "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                  "RGI-CARD", "ABRicate-CARD",
                  "AMRFinderPlus-NCBI", "ABRicate-NCBI",
                  "ResFinder", "ABRicate-ResFinder")

#Same colour but adding texture to the plots
tools_texture <- c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD", 
                   "ABRicate-NCBI", "ABRicate-ResFinder")

# Environments we are not interested in
not_env <- c("amplicon", "isolate", "built-environment")
EN2 <- EN[!EN %in% not_env]
h2 <- c("humans", "mammals", "wastewater", "freshwater", "soil", "marine")  

# Boxplot calculation function
calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  stats <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  iqr <- diff(stats[c(1, 3)])
  
  list(
    ymin = max(min(x, na.rm = TRUE), stats[1] - coef * iqr),
    lower = stats[1],
    middle = stats[2],
    upper = stats[3],
    ymax = min(max(x, na.rm = TRUE), stats[3] + coef * iqr)
  )
}

  
# Define UI for the argCompare application
dashboardPage(
  skin = "green",
  
  dashboardHeader(title = "ARG Detection Tools Comparison"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
      menuItem("Abundance & Diversity", tabName = "abundance", icon = icon("chart-bar")),
      menuItem("Pan & Core Resistome", tabName = "pancore", icon = icon("layer-group")),
      menuItem("Tool Overlap", tabName = "overlap", icon = icon("circle-notch")),
      menuItem("Gene Class Analysis", tabName = "geneclass", icon = icon("dna")),
      menuItem("Data Summary", tabName = "summary", icon = icon("table"))
    ),
    
    hr(),
    
    h4("Filters", style = "padding-left: 15px;"),
    
    pickerInput(
      inputId = "select_tools",
      label = "Select Tools:",
      choices = tools_levels,
      selected = tools_levels,
      multiple = TRUE,
      options = list(`actions-box` = TRUE)
    ),
    
    pickerInput(
      inputId = "select_habitats",
      label = "Select Habitats:",
      choices = EN2,
      selected = c("human gut", "pig gut", "soil"),
      multiple = TRUE,
      options = list(`actions-box` = TRUE)
    ),
    
    sliderInput(
      inputId = "core_threshold",
      label = "Core Resistome Threshold (prevalence):",
      min = 0.1,
      max = 1.0,
      value = 0.5,
      step = 0.1
    ),
    
    numericInput(
      inputId = "min_samples",
      label = "Minimum Samples:",
      value = 450,
      min = 100,
      max = 1000,
      step = 50
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #ecf0f5; }
        .box { margin-bottom: 20px; }
        .info-box { min-height: 90px; }
        .nav-tabs-custom { margin-bottom: 20px; }
      "))
    ),
    
    tabItems(
      
      # Overview Tab
      
      tabItem(
        tabName = "overview",
        fluidRow(
          box(
            title = "About This App",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            HTML("
              <h4>Antibiotic Resistance Gene (ARG) Detection Tool Comparison</h4>
              <p>This interactive application allows you to explore and compare different antibiotic resistance genes detection tools
              across various habitats and environments.</p>
              <h5>Features:</h5>
              <ul>
                <li><strong>Abundance & Diversity:</strong> Compare ARG abundance and diversity metrics across tools and habitats</li>
                <li><strong>Pan & Core Resistome:</strong> Explore pan-resistome and core-resistome patterns</li>
                <li><strong>Tool Overlap:</strong> Visualize agreement and overlap between different detection tools</li>
                <li><strong>Gene Class Analysis:</strong> Examine specific antibiotic resistance gene classes</li>
              </ul>
              <h5>Tools Included:</h5>
              <p>DeepARG, fARGene, ABRicate (with ARG-ANNOT, CARD, MEGARes, NCBI & ResFinder databases), RGI-DIAMOND, AMRFinderPlus, and ResFinder</p>
              <h5>Getting Started:</h5>
              <p>Use the sidebar filters to select tools, habitats, and adjust parameters. Navigate through the tabs to explore different analyses.</p>
            ")
          )
        ),
        
        fluidRow(
          valueBoxOutput("total_samples_box", width = 3),
          valueBoxOutput("total_tools_box", width = 3),
          valueBoxOutput("total_habitats_box", width = 3),
          valueBoxOutput("total_genes_box", width = 3)
        )
      ),
      
      # Abundance & Diversity Tab
      tabItem(
        tabName = "abundance",
        fluidRow(
          box(
            title = "Abundance by Tool and Habitat",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotOutput("plot_abundance", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Diversity (Alpha Diversity)",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotOutput("plot_diversity", height = "600px")
          )
        )
      ),
      
      # Pan & Core Resistome Tab
      tabItem(
        tabName = "pancore",
        fluidRow(
          box(
            title = "Pan vs Core Resistome Size",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotOutput("plot_pancore", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Core Resistome",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotOutput("plot_core_prevalence", height = "400px")
          ),
          box(
            title = "Pan Resistome",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            plotOutput("plot_pan_growth", height = "400px")
          )
        )
      ),
      
      # Tool Overlap Tab
      tabItem(
        tabName = "overlap",
        fluidRow(
          box(
            title = "Tool Agreement Heatmap",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput(
              inputId = "overlap_habitat",
              label = "Select Habitat:",
              choices = EN2,
              selected = "human gut"
            ),
            plotOutput("plot_overlap_heatmap", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Venn Diagram (Select 2-5 tools)",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            pickerInput(
              inputId = "venn_tools",
              label = "Select Tools for Venn:",
              choices = tools_levels,
              selected = tools_levels[1:3],
              multiple = TRUE,
              options = list(`max-options` = 5)
            ),
            plotOutput("plot_venn", height = "500px")
          ),
          box(
            title = "Overlap Statistics",
            width = 6,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            reactableOutput("table_overlap_stats")
          )
        )
      ),
      
      # Gene Class Analysis Tab
      tabItem(
        tabName = "geneclass",
        fluidRow(
          box(
            title = "Gene Class Distribution",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            selectInput(
              inputId = "geneclass_habitat",
              label = "Select Habitat:",
              choices = EN2,
              selected = "human gut"
            ),
            plotOutput("plot_geneclass", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Top Gene Classes",
            width = 12,
            status = "info",
            solidHeader = TRUE,
            collapsible = TRUE,
            reactableOutput("table_geneclass")
          )
        )
      ),
      
      # Data Summary Tab
      tabItem(
        tabName = "summary",
        fluidRow(
          box(
            title = "Sample Summary by Habitat",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            reactableOutput("table_sample_summary")
          )
        ),
        fluidRow(
          box(
            title = "Tool Detection Summary",
            width = 12,
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            reactableOutput("table_tool_summary")
          )
        )
      )
    )
  )
)






