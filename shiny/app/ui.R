library(shiny)
#library(shinydashboard)
library(bslib)
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

source("../../code_R_analysis/helper.R")

# Defining constants
general_size <- 10

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
page_navbar(
  theme = "yeti",
  #title = "How ARG Detection Tools Shape Our View of the Resistome",
  #bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "Introduction", p(intro_page)),
  nav_panel(title = "Fig 1", p(psb1)),
  nav_panel(title = "Fig 2", p(psb1_2)),
  nav_panel(title = "Fig 3", p(psb2)),
  nav_panel(title = "Fig 4", p(psb3)),
  nav_panel(title = "Fig 5", p(psb4)),
  nav_panel(title = "Fig 6", p(psb5)),
  nav_panel(title = "Table S1", p(tab1)),
  nav_panel(title = "Table S2", p(tab2)),
  nav_panel(title = "Table S3", p(tab3)),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("Big Data Biology Lab", href = "https://www.big-data-biology.org")),
    nav_item(tags$a("CMR", href = "https://research.qut.edu.au/cmr/"))
  )
)



# Define server logic 
server <- function(input, output) {
  
  ################################################################################################
  # Table 1
  
  #output$table1 <-  renderTable(df2 %>% arrange(desc(new_level)), striped = TRUE)
  
  # Table 
  
  output$tables1 <- renderReactable({
    reactable(
      t1,
      showPageSizeOptions = TRUE,           # show the dropdown
      pageSizeOptions = c(25, 50, 100, nrow(t1)),  # options in the dropdown
      defaultPageSize = 25
    )
  })
  
  
  output$tables2 <- renderReactable({
    reactable(
      t2,
      showPageSizeOptions = F,           # show the dropdown
      #pageSizeOptions = c(25, 50, 100, nrow(t1)),  # options in the dropdown
      defaultPageSize = nrow(t2)
    )
  })
  
  
  output$tables3 <- renderReactable({
    reactable(
      df2,
      defaultSorted = list(Term_ID = "asc", new_level = "asc"),
      showPageSizeOptions = TRUE,           # show the dropdown
      pageSizeOptions = c(25, 50, 100, nrow(df2)),  # options in the dropdown
      defaultPageSize = 25
    )
  })
  
  # Fig 1
  
  output$unigenes_tool <- renderPlot({
    req(input$tools_unigenes_fig1)
    plot_count_genes_tool(unigenes , input$tools_unigenes_fig1, general_size, pal_10_q, tool_label, tools_levels_2) 
  })
  
  
  output$overlap_heatmap <- renderPlot({ 
    return_heatmap_overalp(overlap_tools, input$tools_unigenes_fig1, general_size, tool_label, tools_levels_2)
  })
  
  # Fig 2
  
  abundance_medians_h2_fixed <- reactive({fix_abundance_diversity_medians_h2(abundance_medians_h2,
                                                                             input$environments_plot, 
                                                                             hab2 = c("Humans","Mammals","Wastewater", "Freshwater","Soil", "Marine"))
  })
  
  diversity_medians_h2_fixed <- reactive({fix_abundance_diversity_medians_h2(diversity_medians_h2,
                                                                             input$environments_plot, 
                                                                             hab2 = c("Humans","Mammals","Wastewater", "Freshwater","Soil", "Marine"))
  })
  
  abundance_medians_h2_fixed_v2 <- reactive({fix_abundance_diversity_medians_version2(abundance_medians_h2,
                                                                                      input$tools_unigenes, 
                                                                                      input$environments_plot)
  })
  
  
  diversity_medians_h2_fixed_v2 <- reactive({fix_abundance_diversity_medians_version2(diversity_medians_h2,
                                                                                      input$tools_unigenes, 
                                                                                      input$environments_plot)
  })
  
  
  output$abundance_tool <- renderPlot({
    plot_total_abundance_diversity_tool(abundance_medians_h2_fixed(), input$tools_unigenes, input$environments_plot , general_size, "abundance", pal_10_q, tool_label, tools_levels_2)
  })
  
  output$diversity_tool <- renderPlot({
    plot_total_abundance_diversity_tool(diversity_medians_h2_fixed(), input$tools_unigenes, input$environments_plot , general_size, "diversity", pal_10_q, tool_label, tools_levels_2)
  })
  
  output$abundance_tool_2 <- renderPlot({
    plot_total_abundance_diversity_tool_version2(abundance_medians_h2_fixed_v2(), input$tools_unigenes, input$environments_plot , general_size, "abundance", pal_6, tool_label, h2, tools_levels_2)
  })
  
  output$diversity_tool_2 <- renderPlot({
    plot_total_abundance_diversity_tool_version2(diversity_medians_h2_fixed_v2(), input$tools_unigenes, input$environments_plot , general_size, "diversity", pal_6, tool_label, h2, tools_levels_2)
  })
  
  ################################################################################################
  # Fig 2
  
  sumcore <- reactive({sum_core_adjust(core, input$cnt_subset, input$threshold_samples )
  })
  
  output$pan_core_tool <- renderPlot({
    pan_resistome_plot(sumpan2, input$tools_levels_fig2, input$environments_plot_fig2, general_size, pal_10_q, tool_label, h2, tools_levels_2)
  })
  
  output$pan_core_tool <- renderPlot({
    plot_fig2(pan_resistome_plot(sumpan2, input$tools_levels_fig2, input$environments_plot_fig2, general_size, pal_10_q, tool_label, h2, tools_levels_2),
              core_resistome_plot(sumcore(), input$tools_levels_fig2, input$environments_plot_fig2, general_size, pal_10_q, tool_label, h2, tools_levels_2))
  })
  
  ################################################################################################   
  # Fig 3
  abundance_env_class <- reactive({abundance_medians_env_class(abundance_class, input$environments_plot_fig3, input$tools_levels_fig3, input$classes_plot_fig3)
  })  
  
  
  unigenes_class <- reactive({get_unigenes_class(unigenes, input$tools_levels_fig3,  input$classes_plot_fig3)
  })
  
  
  output$fig3A <- renderPlot({
    plot_unigenes_env_tool_class(unigenes_class(),  general_size, pal_10_q, tool_label, tools_levels_2, input$tools_levels_fig3)
  })
  
  output$fig3B <- renderPlot({
    plot_abundance_env_tool_class(abundance_env_class(), input$tools_levels_fig3, input$environments_plot_fig3, input$classes_plot_fig3, general_size, pal_10_q, tool_label, tools_levels_2)
  })
  
  output$fig3C <- renderPlot({
    plot_core_env_tool_class(sumcore(), input$tools_levels_fig3, input$environments_plot_fig3, input$classes_plot_fig3, general_size, pal_10_q, tool_label, tools_levels_2)
  })
  
  ################################################################################################
  # Fig 4
  
  
  output$fig4A <- renderPlot({plot_recall_fnr(recall_fnr, input$tool_fig4, input$class_fig4, tool_label, pal_10_q, general_size, tools_levels_2)
  })
  
  output$fig4B <- renderPlot({plot_recall_detailed(recall_fnr, input$tool_fig4, input$class_fig4, tool_label, pal_10_q, general_size, tools_levels_2)
  })
  
  output$fig4C <- renderPlot({plot_fnr_detailed(recall_fnr, input$tool_fig4, input$class_fig4, tool_label, pal_10_q, general_size, tools_levels_2)
  })
  
  ################################################################################################
  # Fig 5
  
  output$fig5A <- renderPlot({
    plot_id_levels(unigenes, tools_levels_2, pal_10_q, tool_label, general_size)
  })  
  
  output$fig5B <- renderPlot({
    merge_recall_fnr_id(unigenes, recall_fnr, input$tool_2_fig5, input$tool_3_fig5, pal_10_q, general_size, tool_label, tools_levels_2)
  })
  
}