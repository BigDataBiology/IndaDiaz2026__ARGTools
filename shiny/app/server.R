server <- function(input, output, session) {
  
  output$plot_count_genes_tool <- renderPlot({
    # req() stops the plot from trying to render if no tools are selected
    req(input$tools_unigenes)
    
    plot_count_genes_tool(
      unigenes = data_list$unigenes, 
      tools_for_figure = input$tools_unigenes, 
      general_size = general_size, 
      pal_10_q = pal_10_q, 
      tool_label = tools_labels , 
      tools_levels = tools_levels, 
      texture = tools_texture) + # which tools have stripes
      theme(panel.background = element_rect(colour = "black", fill = NA)) + 
      theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"))
  })
  
  
  output$plot_alluvial_classes <- renderPlot({
    # req() stops the plot from trying to render if no tools are selected
    req(input$tools_unigenes)
    
    plot_alluvial_classes(data_list$unigenes, 
                          data_list$levels_unigenes, 
                          0.99, 0.005, 
                          input$tools_unigenes, 
                          tools_labels, 
                          tools_levels, 
                          pal_10_q, 
                          general_size, 
                          gene_classes_list) +
      theme(panel.background = element_rect(colour = "black", fill = NA)) 
    
  })
  
}