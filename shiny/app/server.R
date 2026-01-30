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
  
  
  # Core Resistome plot
  pan_core <- reactive({data_list$sumpan2 %>% 
      left_join((sum_core_adjust(data_list$core, input$threshold_samples, input$threshold_proportion) %>% 
                   ungroup() %>% 
                   group_by(tool, habitat) %>% 
                   summarise(core = sum(unigenes))), by = c("tool", "habitat")) %>%
      mutate(core = ifelse(is.na(core), 0, core)) %>% 
      mutate(prop = core / md) %>%
      ungroup() %>% group_by(tool, habitat) %>% 
      mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
      filter(tool %in% input$tool_pan_core,
             habitat %in% input$environment_pan_core) %>%
      mutate(tool = factor(as.character(tool),
                 levels = tools_levels[tools_levels %in% input$tool_pan_core]))})
  
  output$plot_pan_core_resistome <- renderPlot({
    tools_order <- match(input$tool_pan_core, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    pal_figure <- pal_10_q[tools_order]
    tools_labels_figure <- tools_labels[tools_order]
    
    
    pan_core() %>% select(!c(md, sd)) %>% pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
      mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
      mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
      mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>% 
      ggplot(aes(x = habitat, y =  value)) +
        geom_jitter(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 2.5, width = 0.7, height = 0) + 
        facet_grid(metric ~ habitat, scales = "free") +
        scale_fill_manual(values = pal_figure, labels = lab_fn(tools_labels_figure), name = NULL) +
        scale_shape_manual(values = c(21, 24)) +
        guides(
          fill = guide_legend(
            override.aes = list(
              shape = shape_tools,
              fill  = pal_figure)), shape = "none") +
        theme_minimal() +
        xlab("") +
        ylab("ARGs") + 
        theme(
          legend.position = "bottom",
          strip.text.x   = element_text(size = general_size),
          #strip.text.y   = element_blank(),
          legend.text = element_text(size = general_size),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          panel.spacing = unit(0, "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = general_size),
          panel.border = element_blank(),   
          panel.background = element_rect(colour = "black", fill = NA))  
  })
  
}