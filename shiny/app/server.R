server <- function(input, output, session) {
  
  input_data_unigenes <- reactive({data_list$unigenes %>% 
    filter(!(tool %in% c("DeepARG", "RGI-DIAMOND") &  id < input$threshold_unigenes_id))})
  
  output$plot_count_genes_tool <- renderPlot({
    
    input_data <- data_list$unigenes %>% 
      filter(!(tool %in% c("DeepARG", "RGI-DIAMOND") &  id < input$threshold_unigenes_id))
    
    plot_count_genes_tool(
      unigenes = input_data_unigenes(), 
      tools_for_figure = input$tools_unigenes, 
      general_size = general_size, 
      pal_10_q = pal_10_q, 
      tool_label = tools_labels , 
      tools_levels = tools_levels, 
      texture = tools_texture) + 
      theme(panel.background = element_rect(colour = "black", fill = NA)) + 
      theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"))
  })
  
  output$plot_alluvial_classes <- renderPlot({

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
          bind_rows(data_list$sumpan2_60) %>%  
          left_join(
            (
              sum_core_adjust(
                (data_list$core %>% 
                  filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                  bind_rows(data_list$core60)), 
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
          bind_rows(data_list$sumpan2_70) %>%  
          left_join(
            (
              sum_core_adjust(
                (data_list$core %>% 
                  filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                  bind_rows(data_list$core70)), 
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
          bind_rows(data_list$sumpan2_80) %>%  
          left_join(
            (
              sum_core_adjust(
                (data_list$core %>% 
                  filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
                  bind_rows(data_list$core80)), 
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
        data_list$sumpan2  %>%  
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
  
  
  ### ABUNDANCE
  
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
      pal_10_q = pal_10_q, # pallet
      metric = "abundance", # metric (abundance or diversity)
      sd = 2025, # seed to plot random samples in the distribution 
      obs = 200,  # number of samples to plot as dots per environment
      texture = tools_texture, # texture for repeated color 
      tools_levels = tools_levels) + 
      theme(legend.position = "none")
    
  })
  
  output$plot_diversity <- renderPlot({
    
    plot_total_abundance_diversity_new_version(
      dataset = abundance_tool_sample_reactive(), # 
      tools_labels = tools_labels,  #
      tools_to_plot = input$tool_abundance,  #
      environments_plot = input$environment_abundance, # habitats to plot (aggregated humans and mammals)
      general_size = general_size, # font size
      pal_10_q = pal_10_q, # pallet
      metric = "diversity", # metric (abundance or diversity)
      sd = 2025, # seed to plot random samples in the distribution 
      obs = 200,  # number of samples to plot as dots per environment
      texture = tools_texture, # texture for repeated color 
      tools_levels = tools_levels) + 
      theme(legend.position = "none")
    
  })

  
  abundance_class_reactice <- reactive({ 
    if(input$threshold_abundance_id == 60.0) {
      data_list$abundance_class %>% 
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance_class60) %>%
        filter(tool %in% input$tool_abundance)
      
    } else if(input$threshold_abundance_id == 70.0){
      data_list$abundance_class %>% 
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance_class70) %>%
        filter(tool %in% input$tool_abundance)
      
    } else if(input$threshold_abundance_id == 80.0) {
      data_list$abundance_class %>% 
        filter(!tool %in% c("DeepARG","RGI-DIAMOND")) %>% 
        bind_rows(data_list$abundance_class80) %>%
        filter(tool %in% input$tool_abundance)
      
    } else {data_list$abundance_class %>%
        filter(tool %in% input$tool_abundance)}
    
  })
  
  output$plot_abundance_class  <- renderPlot({
    plot_abundance_class_more_environments(abundance_class_reactice(), 
          input$environment_abundance, general_size , 
          pal_10_q, 
          input$abundance_genes, 
          data_type = "abundance", 
          other = input$plot_other,
          tools_levels, 
          input$tool_abundance, 
          tools_texture,
          pattern_density = pattern_density,
          pattern_spacing = pattern_spacing,
          pattern_fill = pattern_fill,
          pattern_size = pattern_size) +
      theme(legend.position = "none")
  })
  
  output$plot_diversity_class  <- renderPlot({
    plot_abundance_class_more_environments(abundance_class_reactice(), 
      input$environment_abundance, general_size , 
      pal_10_q, 
      input$abundance_genes, 
      data_type = "diversity", 
      other = input$plot_other,
      tools_levels, 
      input$tool_abundance, 
      tools_texture,
      pattern_density = pattern_density,
      pattern_spacing = pattern_spacing,
      pattern_fill = pattern_fill,
      pattern_size = pattern_size) +
      theme(legend.position = "none")
  })
  
  output$plot_diversity_class_legend  <- renderPlot({
    grid.draw(g_legend(plot_abundance_class_more_environments(abundance_class_reactice(), 
        input$environment_abundance, general_size , 
        pal_10_q, 
        input$abundance_genes, 
        data_type = "diversity", 
        other = input$plot_other,
        tools_levels, 
        input$tool_abundance, 
        tools_texture,
        pattern_density = 0.01,
        pattern_spacing = 0.005,
        pattern_fill = "white",
        pattern_size = 0.4) +
      theme(legend.position = "bottom")))
  })
  
  
  ### OVERLAPS 

  recall_reactive <- reactive({
    if(input$threshold_overlap_id == 60.0){
      data_list$recall_fnr60  %>% 
        mutate(new_level = as.character(new_level)) %>% 
        filter(new_level %in% input$overlap_genes) %>% 
        filter(tool_ref %in% input$tool_overlap,
               tool_comp %in% input$tool_overlap_calc) %>%
        mutate(tool_ref = factor(as.character(tool_ref), 
                                 levels = tools_levels[tools_levels %in% 
                                                         input$tool_overlap])) %>% 
        mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
        mutate(new_level = factor(new_level, levels = input$overlap_genes))
      
    } else if(input$threshold_overlap_id == 70.0){
      data_list$recall_fnr70  %>% 
        mutate(new_level = as.character(new_level)) %>% 
        filter(new_level %in% input$overlap_genes) %>% 
        filter(tool_ref %in% input$tool_overlap,
               tool_comp %in% input$tool_overlap_calc) %>%
        mutate(tool_ref = factor(as.character(tool_ref), 
                                 levels = tools_levels[tools_levels %in% 
                                                         input$tool_overlap])) %>% 
        mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
        mutate(new_level = factor(new_level, levels = input$overlap_genes))
      
    } else if(input$threshold_overlap_id == 80.0){
      data_list$recall_fnr80  %>% 
        mutate(new_level = as.character(new_level)) %>% 
        filter(new_level %in% input$overlap_genes) %>% 
        filter(tool_ref %in% input$tool_overlap,
               tool_comp %in% input$tool_overlap_calc) %>%
        mutate(tool_ref = factor(as.character(tool_ref), 
                                 levels = tools_levels[tools_levels %in% 
                                                         input$tool_overlap])) %>% 
        mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
        mutate(new_level = factor(new_level, levels = input$overlap_genes))

    } else {
      data_list$recall_fnr  %>% 
        mutate(new_level = as.character(new_level)) %>% 
        filter(new_level %in% input$overlap_genes) %>% 
        filter(tool_ref %in% input$tool_overlap,
               tool_comp %in% input$tool_overlap_calc) %>%
        mutate(tool_ref = factor(as.character(tool_ref), 
                                 levels = tools_levels[tools_levels %in% 
                                                         input$tool_overlap])) %>% 
        mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
        mutate(new_level = factor(new_level, levels = input$overlap_genes))
    }
  })
    
  output$overlap_cstc  <- renderPlot({

    recall_reactive() %>%
      ggplot(aes(x = new_level, fill = tool_ref, y = recall)) +
        geom_violin() +
        geom_jitter(aes(color = tool_comp, shape = texture), #color = "black", 
                    stroke = 1, size = 2.5, width = 0.1,  height = 0) + 
        scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
        scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                          labels = tools_levels[tools_levels %in% input$tool_overlap]) +
        scale_color_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
                           labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
        scale_x_discrete(labels = function(x) {
          x <- gsub("-", "-\n", x)
          x <- gsub(" ", "\n", x)
          x}) + 
        scale_shape_manual(values = c("no" = 16, "yes" = 17)) +
        facet_grid(tool_ref ~ new_level, scales = "free_x") +
        scale_y_continuous(limits = c(-0.1,1.1), 
                           breaks = seq(0, 1, length.out = 3),
                           labels = scales::label_number()) + 
        ylab("CSTC") + 
        xlab("") + 
        theme(
          legend.position = "none",
          legend.text = element_text(size = general_size ),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          #panel.spacing = unit(0, "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_text(size = general_size),
          panel.background = element_rect(colour = "black", fill = NA),
          strip.text = element_blank(),
          strip.background = element_blank())

  })
  
  
  output$overlap_cstc_summary  <- renderPlot({
    #pal2 <- pal_10_q[]
    #pal3 <- pal_10_q[input$tool_overlap_calc %in% tool_levels]
    
    recall_reactive() %>%
      filter(!is.na(recall)) %>% 
      ungroup() %>% 
      group_by(tool_ref, new_level, texture) %>% 
      summarise(recall = median(recall)) %>%
      ggplot(aes(x = "Class medians", fill = tool_ref, y = recall)) +
        geom_violin() +
        #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
        geom_jitter(aes(shape = texture), size = 1.5, width = 0.4, height = 0) + 
        scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
        scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                          labels = tools_levels[tools_levels %in% input$tool_overlap]) +
        #scale_color_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
        #                   labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
        facet_grid(tool_ref ~ ., scales = "free_x") +
        scale_y_continuous(limits = c(-0.2,1.2), 
                           breaks = seq(0, 1, length.out = 3),
                           labels = scales::label_number()) +
        scale_x_discrete(labels = function(x) {
          x <- gsub("-", "-\n", x)
          x <- gsub(" ", "\n", x)
          x}) + 
        ylab("") + 
        xlab("") + 
        theme(
          legend.position = "none",
          legend.text = element_text(size = general_size ),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          panel.spacing = unit(0, "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_blank(),
          panel.background = element_rect(colour = "black", fill = NA),
          strip.text = element_text(size = general_size, face = "bold"),
          strip.background = element_blank())
    
  })
  
  
  output$overlap_csno  <- renderPlot({
    
    recall_reactive() %>%
      ggplot(aes(x = new_level, fill = tool_ref, y = fnr)) +
        geom_violin() +
        geom_jitter(aes(color = tool_comp, shape = texture), #color = "black", 
                    stroke = 1, size = 2.5, width = 0.1,  height = 0) + 
        scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
        scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                          labels = tools_levels[tools_levels %in% input$tool_overlap]) +
        scale_color_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
                           labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
        scale_x_discrete(labels = function(x) {
          x <- gsub("-", "-\n", x)
          x <- gsub(" ", "\n", x)
          x}) + 
        scale_shape_manual(values = c("no" = 16, "yes" = 17)) +
        facet_grid(tool_ref ~ new_level, scales = "free_x") +
        scale_y_continuous(limits = c(-0.1,1.1), 
                           breaks = seq(0, 1, length.out = 3),
                           labels = scales::label_number()) + 
        ylab("CSNO") + 
        xlab("") + 
        theme(
          legend.position = "none",
          legend.text = element_text(size = general_size ),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          #panel.spacing = unit(0, "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_text(size = general_size),
          panel.background = element_rect(colour = "black", fill = NA),
          strip.text = element_blank(),
          strip.background = element_blank())
      
  })
  
  
  output$overlap_csno_summary  <- renderPlot({
    #pal2 <- pal_10_q[]
    #pal3 <- pal_10_q[input$tool_overlap_calc %in% tool_levels]
    
    recall_reactive() %>%
      filter(!is.na(recall)) %>% 
      ungroup() %>% 
      group_by(tool_ref, new_level, texture) %>% 
      summarise(fnr = median(fnr)) %>%
      ggplot(aes(x = "Class medians", fill = tool_ref, y = fnr)) +
        geom_violin() +
        #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
        geom_jitter(aes(shape = texture), size = 1.5, width = 0.4, height = 0) + 
        scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
        scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                          labels = tools_levels[tools_levels %in% input$tool_overlap]) +
        #scale_color_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
        #                   labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
        facet_grid(tool_ref ~ ., scales = "free_x") +
        scale_y_continuous(limits = c(-0.2,1.2), 
                           breaks = seq(0, 1, length.out = 3),
                           labels = scales::label_number()) +
        scale_x_discrete(labels = function(x) {
          x <- gsub("-", "-\n", x)
          x <- gsub(" ", "\n", x)
          x}) + 
        ylab("") + 
        xlab("") + 
        theme(
          legend.position = "none",
          legend.text = element_text(size = general_size ),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          panel.spacing = unit(0, "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_blank(),
          panel.background = element_rect(colour = "black", fill = NA),
          strip.text = element_text(size = general_size, face = "bold"),
          strip.background = element_blank())
      
  })
  
  
  output$overlap_legend <- renderPlot({
    
    tools_order <- match(input$tool_overlap_calc, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]

    grid.draw(g_legend(recall_reactive() %>%
      ggplot(aes(x = recall + fnr, y = recall + fnr)) +
      geom_jitter(aes(shape = texture, fill = tool_comp), #color = "black", 
                  size = 2.5) + 
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap_calc], 
                         labels = tools_levels[tools_levels %in% input$tool_overlap_calc]) +
      scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
      guides(
        fill = guide_legend(
          override.aes = list(
            shape = shape_tools,
            #color = NA,   
            #stroke = 0,
            fill  = pal_10_q[tools_levels %in% input$tool_overlap_calc])), shape = "none") +
      labs(fill = "") +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = general_size ),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        #panel.spacing = unit(0, "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        panel.background = element_rect(colour = "black", fill = NA),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(), 
        legend.background = element_blank())))
  })
}


