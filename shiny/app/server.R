server <- function(input, output, session) {
  
  input_data_unigenes <- reactive({data_list$unigenes %>% 
    filter(!(tool %in% c("DeepARG", "RGI-DIAMOND") &  id < input$threshold_unigenes_id))})
  
  output$plot_count_genes_tool <- renderPlot({
    # req() stops the plot from trying to render if no tools are selected
    req(input$tools_unigenes)
    
    input_data <- data_list$unigenes %>% 
      filter(!(tool %in% c("DeepARG", "RGI-DIAMOND") &  id < input$threshold_unigenes_id))
    
    plot_count_genes_tool(
      unigenes = input_data_unigenes(), 
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

    plot_alluvial_classes(unigenes = input_data_unigenes(), 
                          levels_unigenes = data_list$levels_unigenes, 
                          threshold_plot = 0.99, remove_class_threshold = 0.005, 
                          tools_to_plot = input$tools_unigenes, 
                          tools_labels = tools_labels, 
                          tools_factors = tools_levels, 
                          pal_10_q = pal_10_q, 
                          general_size = general_size, 
                          gene_classes_list = gene_classes_list) +
      theme(panel.background = element_rect(colour = "black", fill = NA)) 
    
  })
  

  # Core Resistome plot
  pan_core <- reactive({ 
      
      if(input$threshold_pan_core_id == 60.0) {
        data_list$sumpan2 %>% 
          filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
          bind_rows(data_list$sumpan2_60) %>%  ## here
          left_join(
            (
              sum_core_adjust(
                (data_list$core %>% 
                  filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                  bind_rows(data_list$core60)), ## here
                input$threshold_samples, input$threshold_proportion) %>% 
                       ungroup() %>% 
                       group_by(tool, habitat) %>% 
                       summarise(core = sum(unigenes))), 
            by = c("tool", "habitat")) %>%
          mutate(core = ifelse(is.na(core), 0, core)) %>% 
          mutate(prop = core / md) %>%
          ungroup() %>% group_by(tool, habitat) %>% 
          mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
          filter(tool %in% input$tool_pan_core,
                 habitat %in% input$environment_pan_core) %>%
          mutate(tool = factor(as.character(tool),
                               levels = tools_levels[tools_levels %in% input$tool_pan_core]))
        
      } else if(input$threshold_pan_core_id == 70.0) {
        data_list$sumpan2 %>% 
          filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
          bind_rows(data_list$sumpan2_70) %>%  ## here
          left_join(
            (
              sum_core_adjust(
                (data_list$core %>% 
                  filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                  bind_rows(data_list$core70)), ## here
                input$threshold_samples, input$threshold_proportion) %>% 
                ungroup() %>% 
                group_by(tool, habitat) %>% 
                summarise(core = sum(unigenes))), 
            by = c("tool", "habitat")) %>%
          mutate(core = ifelse(is.na(core), 0, core)) %>% 
          mutate(prop = core / md) %>%
          ungroup() %>% group_by(tool, habitat) %>% 
          mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
          filter(tool %in% input$tool_pan_core,
                 habitat %in% input$environment_pan_core) %>%
          mutate(tool = factor(as.character(tool),
                               levels = tools_levels[tools_levels %in% input$tool_pan_core]))
        
      } else if(input$threshold_pan_core_id == 80.0) {
        data_list$sumpan2 %>% 
          filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
          bind_rows(data_list$sumpan2_80) %>%  ## here
          left_join(
            (
              sum_core_adjust(
                (data_list$core %>% 
                  filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                  bind_rows(data_list$core80)), ## here
                input$threshold_samples, input$threshold_proportion) %>% 
                ungroup() %>% 
                group_by(tool, habitat) %>% 
                summarise(core = sum(unigenes))), 
            by = c("tool", "habitat")) %>%
          mutate(core = ifelse(is.na(core), 0, core)) %>% 
          mutate(prop = core / md) %>%
          ungroup() %>% group_by(tool, habitat) %>% 
          mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
          filter(tool %in% input$tool_pan_core,
                 habitat %in% input$environment_pan_core) %>%
          mutate(tool = factor(as.character(tool),
                               levels = tools_levels[tools_levels %in% input$tool_pan_core]))
        
        
      } else {
        data_list$sumpan2  %>%  ## here
          left_join(
            (
              sum_core_adjust(
                data_list$core,
                input$threshold_samples, 
                input$threshold_proportion) %>% 
              ungroup() %>% 
              group_by(tool, habitat) %>% 
              summarise(core = sum(unigenes))), 
            by = c("tool", "habitat")) %>%
          mutate(core = ifelse(is.na(core), 0, core)) %>% 
          mutate(prop = core / md) %>%
          ungroup() %>% group_by(tool, habitat) %>% 
          mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
          filter(tool %in% input$tool_pan_core,
                 habitat %in% input$environment_pan_core) %>%
          mutate(tool = factor(as.character(tool),
                               levels = tools_levels[tools_levels %in% input$tool_pan_core]))
      }
})

    
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
        geom_jitter(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 2.5, width = 0.5, height = 0) + 
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
  
  
  abundance_tool_sample_reactive <- reactive({ 
    if(input$threshold_abundance_id == 60.0) {
      data_list$abundance %>%
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance60) %>%
        group_by(tool, sample, habitat, habitat2) %>%  
        summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
        ungroup() %>% 
        complete(sample, tool) %>% # complete with NAs
        
        left_join((data_list$abundance %>%
                     filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                     bind_rows(data_list$abundance60)) %>% select(sample, habitat, habitat2) %>% 
                    distinct(), by = "sample") %>% # get habitat and habitat2 
        mutate(habitat  = coalesce(habitat.x, habitat.y), 
               habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
        select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
        mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
        mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
        arrange(tool, sample)
      
    } else if(input$threshold_abundance_id == 70.0) {
      
      data_list$abundance %>%
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance70) %>%
        group_by(tool, sample, habitat, habitat2) %>%  
        summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
        ungroup() %>% 
        complete(sample, tool) %>% # complete with NAs
        
        left_join((data_list$abundance %>%
                     filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                     bind_rows(data_list$abundance70)) %>% select(sample, habitat, habitat2) %>% 
                    distinct(), by = "sample") %>% # get habitat and habitat2 
        mutate(habitat  = coalesce(habitat.x, habitat.y), 
               habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
        select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
        mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
        mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
        arrange(tool, sample)
      
    } else if (input$threshold_abundance_id == 80.0) {
      
      data_list$abundance %>%
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance80) %>%
        group_by(tool, sample, habitat, habitat2) %>%  
        summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
        ungroup() %>% 
        complete(sample, tool) %>% # complete with NAs
        
        left_join((data_list$abundance %>%
                     filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                     bind_rows(data_list$abundance80)) %>% select(sample, habitat, habitat2) %>% 
                    distinct(), by = "sample") %>% # get habitat and habitat2 
        mutate(habitat  = coalesce(habitat.x, habitat.y), 
               habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
        select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
        mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
        mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
        arrange(tool, sample)
      
    } else {
      data_list$abundance %>%
        group_by(tool, sample, habitat, habitat2) %>%  
        summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
        ungroup() %>% 
        complete(sample, tool) %>% # complete with NAs
        left_join(data_list$abundance %>% select(sample, habitat, habitat2) %>% 
                    distinct(), by = "sample") %>% # get habitat and habitat2 
        mutate(habitat  = coalesce(habitat.x, habitat.y), 
               habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
        select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
        mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
        mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
        arrange(tool, sample)
    }
  })
      
  output$plot_abundance <- renderPlot({

    plot_total_abundance_diversity_new_version(
      dataset = abundance_tool_sample_reactive(), # 
      tools_labels = tools_labels,  #
      tools_to_plot = input$tool_abundance,  #
      environments_plot = input$environment_abundance, # habitats to plot (aggregated humans and mammals)
      general_size = general_size, # font size
      pal_10_q = pal_10_q , # pallet
      metric = "abundance", # metric (abundance or diversity)
      sd = 2025, # seed to plot random samples in the distribution 
      obs = 200,  # number of samples to plot as dots per environment
      texture = tools_texture) + # texture for repeated color 
      theme(legend.position = "none")
    
  })
  
  output$plot_diversity <- renderPlot({
    
    plot_total_abundance_diversity_new_version(
      dataset = abundance_tool_sample_reactive(), # 
      tools_labels = tools_labels,  #
      tools_to_plot = input$tool_abundance,  #
      environments_plot = input$environment_abundance, # habitats to plot (aggregated humans and mammals)
      general_size = general_size, # font size
      pal_10_q = pal_10_q , # pallet
      metric = "diversity", # metric (abundance or diversity)
      sd = 2025, # seed to plot random samples in the distribution 
      obs = 200,  # number of samples to plot as dots per environment
      texture = tools_texture) + # texture for repeated color 
      theme(legend.position = "none")
    
  })

  
  abundance_class_reactice <- reactive({ 
    if(input$threshold_abundance_id == 60.0) {
      data_list$abundance_class %>% 
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance_class60)
      
    } else if(input$threshold_abundance_id == 70.0){
      data_list$abundance_class %>% 
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance_class70)
      
    } else if(input$threshold_abundance_id == 80.0) {
      data_list$abundance_class %>% 
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance_class80)
      
    } else {data_list$abundance_class}
    
  })
  
  output$plot_abundance_class  <- renderPlot({
    plot_abundance_class_more_environments(abundance_class_reactice(), 
      input$environment_abundance, general_size , pal_10_q, 
      input$abundance_genes, data_type = "abundance", other = input$plot_other) +
      theme(legend.position = "none")
  })
  
  output$plot_diversity_class  <- renderPlot({
    plot_abundance_class_more_environments(abundance_class_reactice(), 
      input$environment_abundance, general_size , pal_10_q, 
      input$abundance_genes, data_type = "diversity", other = input$plot_other) +
      theme(legend.position = "none")
  })
  
}


