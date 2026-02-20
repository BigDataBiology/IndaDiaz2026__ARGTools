server <- function(input, output, session) {


  ## Caching Unigenes for a fast loading
  unigenes_filtered <- reactive({
    data_list$unigenes %>% filter(!(tool %in% c("DeepARG", "RGI-DIAMOND") & id < input$threshold_unigenes_id))
  })
  
  output$plot_count_genes_tool <- renderPlot({
    
    plot_count_genes_tool(
      unigenes = unigenes_filtered(), 
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
    
    plot_alluvial_classes(unigenes = unigenes_filtered(), 
                          levels_unigenes = data_list$levels_unigenes, 
                          threshold_plot = 0.99, remove_class_threshold = 0.005, 
                          tools_to_plot = input$tools_unigenes, 
                          tools_labels = tools_labels, 
                          tools_factors = tools_levels, 
                          pal_10_q = pal_10_q, 
                          general_size = general_size, 
                          gene_classes_list = gene_classes_list
                          ) +
      theme(panel.background = element_rect(colour = "black", fill = NA))
    
  })
  
  # Core Resistome plot
  
  pan_core <- reactive({
    req(input$threshold_pan_core_id, input$threshold_samples, input$threshold_proportion,
        input$tool_pan_core, input$environment_pan_core)
    
    tools_excl <- c("DeepARG", "RGI-DIAMOND")
    thr <- as.character(input$threshold_pan_core_id)
    
    # Pick threshold-specific "extra" tables (or NULL/default)
    extra_sumpan2 <- switch(
      thr,
      "60" = data_list$sumpan2_60,
      "70" = data_list$sumpan2_70,
      "80" = data_list$sumpan2_80,
      NULL
    )
    
    extra_core <- switch(
      thr,
      "60" = data_list$core60,
      "70" = data_list$core70,
      "80" = data_list$core80,
      NULL
    )
    
    # Base tables:
    sumpan2_df <- if (is.null(extra_sumpan2)) {
      data_list$sumpan2
    } else {
      dplyr::bind_rows(
        data_list$sumpan2 %>% dplyr::filter(!tool %in% tools_excl),
        extra_sumpan2
      )
    }
    
    core_df <- if (is.null(extra_core)) {
      data_list$core
    } else {
      dplyr::bind_rows(
        data_list$core %>% dplyr::filter(!tool %in% tools_excl),
        extra_core
      )
    }
    
    # Compute core summary once
    core_sum <- sum_core_adjust(
      core_df,
      input$threshold_samples,
      input$threshold_proportion
    ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(tool, habitat) %>%
      dplyr::summarise(core = sum(unigenes, na.rm = TRUE), .groups = "drop")
    
    # Apply join + common transforms once
    sumpan2_df %>%
      dplyr::left_join(core_sum, by = c("tool", "habitat")) %>%
      dplyr::mutate(
        core = tidyr::replace_na(core, 0),
        prop = core / md,
        texture = ifelse(tool %in% tools_texture, "yes", "no")
      ) %>%
      dplyr::filter(
        tool %in% input$tool_pan_core,
        habitat %in% input$environment_pan_core
      ) %>%
      dplyr::mutate(
        tool = factor(as.character(tool),
                      levels = tools_levels[tools_levels %in% input$tool_pan_core])
      )
  }) %>% bindCache(
    input$threshold_pan_core_id,
    input$threshold_samples,
    input$threshold_proportion,
    sort(input$tool_pan_core),
    sort(input$environment_pan_core)
  )
  
  
  output$plot_pan_core_resistome <- renderPlot({
    req(input$tool_pan_core, input$environment_pan_core)
    
    tools_order <- match(input$tool_pan_core, tools_levels)
    
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    
    pal_figure <- pal_10_q[tools_order]
    tools_labels_figure <- tools_labels[tools_order]
    
    pan_core() %>%
      dplyr::select(!c(md, sd)) %>%
      tidyr::pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>%
      dplyr::mutate(
        metric = dplyr::recode(metric, mn = "Pan-resistome", core = "Core-resistome"),
        metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))
      ) %>%
      ggplot2::ggplot(ggplot2::aes(x = habitat, y = value)) +
      ggplot2::geom_jitter(
        ggplot2::aes(fill = tool, shape = texture),
        color = "black", stroke = 0.3, size = 2.5,
        width = 0.5, height = 0
      ) +
      ggplot2::facet_grid(metric ~ habitat, scales = "free") +
      ggplot2::scale_fill_manual(values = pal_figure, labels = lab_fn(tools_labels_figure), name = NULL) +
      ggplot2::scale_shape_manual(values = c(21, 24)) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          override.aes = list(shape = shape_tools, fill = pal_figure)
        ),
        shape = "none"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::xlab("") +
      ggplot2::ylab("ARGs") +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text.x   = ggplot2::element_text(size = general_size),
        legend.text    = ggplot2::element_text(size = general_size),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = grid::unit(0, "pt"),
        title = ggplot2::element_text(size = general_size + 2, face = "bold"),
        axis.title = ggplot2::element_text(size = general_size + 1, face = "bold"),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = general_size),
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(colour = "black", fill = NA)
      )
  })
  
  
  
  ### ABUNDANCE
  
  abundance_tool_sample_reactive <- reactive({
    req(input$threshold_abundance_id)
    
    switch(
      as.character(input$threshold_abundance_id),
      "60" = abundance_prepped[["60"]],
      "70" = abundance_prepped[["70"]],
      "80" = abundance_prepped[["80"]],
      abundance_prepped[["default"]]
    )
  }) %>% bindCache(input$threshold_abundance_id)
  
  
  # Filter by selected tools + environments
  abundance_filtered <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    abundance_tool_sample_reactive() %>%
      dplyr::filter(
        tool %in% input$tool_abundance,
        habitat %in% input$environment_abundance
      )
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(input$tool_abundance),
    sort(input$environment_abundance)
  )
  
  

  plot_abundance_reactive <- reactive({
    req(input$tool_abundance, input$environment_abundance, input$threshold_abundance_id)
    
    df <- abundance_filtered() %>%
      dplyr::mutate(tool = factor(tool, levels = input$tool_abundance))
    
    plot_total_abundance_diversity_new_version(
      dataset = df,
      tools_labels = tools_labels,
      tools_to_plot = input$tool_abundance,       
      environments_plot = input$environment_abundance,
      general_size = general_size,
      pal_10_q = pal_10_q,
      metric = "abundance",
      sd = 2025,
      obs = 200,
      texture = tools_texture,
      tools_levels = tools_levels
    ) +
      theme(legend.position = "none")
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(input$tool_abundance),                      
    sort(input$environment_abundance)
  )
  

  ## Adding this so it can be rendered separately well in the submenus
  output$plot_abundance <- renderPlot(plot_abundance_reactive())
  output$plot_abundance_overview <- renderPlot(plot_abundance_reactive())
  
  
  plot_diversity_reactive <- reactive({
    req(input$tool_abundance, input$environment_abundance, input$threshold_abundance_id)
    
    df <- abundance_filtered() %>%
      dplyr::mutate(tool = factor(tool, levels = input$tool_abundance))
    
    plot_total_abundance_diversity_new_version(
      dataset = df,
      tools_labels = tools_labels,
      tools_to_plot = input$tool_abundance,         
      environments_plot = input$environment_abundance,
      general_size = general_size,
      pal_10_q = pal_10_q,
      metric = "diversity",
      sd = 2025,
      obs = 200,
      texture = tools_texture,
      tools_levels = tools_levels
    ) +
      theme(legend.position = "none")
  }) %>% bindCache(
    input$threshold_abundance_id,
    sort(input$tool_abundance),                      
    sort(input$environment_abundance)
  )

  ## Adding this so it can be rendered separately well in the submenus
  output$plot_diversity <- renderPlot(plot_diversity_reactive())
  output$plot_diversity_overview <- renderPlot(plot_diversity_reactive())
  
  ##Suspend any hidden plots
  outputOptions(output, "plot_abundance", suspendWhenHidden = TRUE)
  outputOptions(output, "plot_abundance_overview", suspendWhenHidden = TRUE)
  outputOptions(output, "plot_diversity", suspendWhenHidden = TRUE)
  outputOptions(output, "plot_diversity_overview", suspendWhenHidden = TRUE)

  
  # Median abundance and diversity
  abundance_class_reactice <- (
    reactive({
      req(input$threshold_abundance_id, input$tool_abundance)
      
      base <- data_list$abundance_class
      
      extra <- switch(
        as.character(input$threshold_abundance_id),
        "60" = data_list$abundance_class60,
        "70" = data_list$abundance_class70,
        "80" = data_list$abundance_class80,
        NULL
      )
      
      dplyr::bind_rows(base, extra) %>%
        dplyr::filter(tool %in% input$tool_abundance) %>%                        
        dplyr::mutate(tool = factor(tool, levels = input$tool_abundance))       
    })
  ) %>% bindCache(
    input$threshold_abundance_id,
    sort(input$tool_abundance)
  )

  
  plot_abundance_class_reactive <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    sel_tools <- intersect(tools_levels, input$tool_abundance)
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    plot_abundance_class_more_environments(
      abundance_class_reactice(),
      input$environment_abundance, general_size,
      pal_sel,                      
      input$abundance_genes,
      data_type = "abundance",
      other = input$plot_other,
      tools_levels,
      sel_tools,                     
      tools_texture,
      pattern_density = pattern_density,
      pattern_spacing = pattern_spacing,
      pattern_fill = pattern_fill,
      pattern_size = pattern_size
    ) +
      theme(legend.position = "none")
  })
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$plot_abundance_class <- renderPlot(plot_abundance_class_reactive())
  output$plot_abundance_class_overview <- renderPlot(plot_abundance_class_reactive())
  
  
  plot_diversity_class_reactive <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    sel_tools <- intersect(tools_levels, input$tool_abundance)
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    plot_abundance_class_more_environments(
      abundance_class_reactice(),
      input$environment_abundance, general_size,
      pal_sel,                        
      input$abundance_genes,
      data_type = "diversity",
      other = input$plot_other,
      tools_levels,
      sel_tools,
      tools_texture,
      pattern_density = pattern_density,
      pattern_spacing = pattern_spacing,
      pattern_fill = pattern_fill,
      pattern_size = pattern_size
    ) +
      theme(legend.position = "none")
  })
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$plot_diversity_class <- renderPlot(plot_diversity_class_reactive())
  output$plot_diversity_class_overview <- renderPlot(plot_diversity_class_reactive())
  
  plot_diversity_class_legend <- reactive({
    req(input$threshold_abundance_id, input$environment_abundance, input$tool_abundance)

    sel_tools <- intersect(tools_levels, input$tool_abundance)
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)

    p <- plot_abundance_class_more_environments(
      abundance_class_reactice(),
      input$environment_abundance, general_size,
      pal_sel,       #So the input of the user reflects in the plot                 
      input$abundance_genes,
      data_type = "diversity",
      other = input$plot_other,
      tools_levels,
      sel_tools,
      tools_texture,
      pattern_density = 0.01,
      pattern_spacing = 0.005,
      pattern_fill = "white",
      pattern_size = 0.4
    ) +
      theme(legend.position = "bottom")

    g_legend(p)
  })
  
  
  ## Adding this so the legend can be rendered in the plot
  output$plot_diversity_class_legend <- renderPlot({
    grid::grid.newpage()
    grid::grid.draw(plot_diversity_class_legend())
  })
  
  output$plot_diversity_class_legend_overview <- renderPlot({
    grid::grid.newpage()
    grid::grid.draw(plot_diversity_class_legend())
  })
  

  plot_abundance_class_legend <- reactive({
    req(input$threshold_abundance_id,
        input$environment_abundance,
        input$tool_abundance)
    
    sel_tools <- intersect(tools_levels, input$tool_abundance)
    pal_sel   <- pal_for_tools(sel_tools, tools_levels, pal_10_q)
    
    p <- plot_abundance_class_more_environments(
      abundance_class_reactice(),
      input$environment_abundance,
      general_size,
      pal_sel,
      input$abundance_genes,
      data_type = "abundance",
      other = input$plot_other,
      tools_levels,
      sel_tools,
      tools_texture,
      pattern_density = 0.01,
      pattern_spacing = 0.005,
      pattern_fill = "white",
      pattern_size = 0.4
    ) +
      theme(legend.position = "bottom")
    
    g_legend(p)
  })
  
  ## Adding this so the legend can be rendered in the plot
  output$plot_abundance_class_legend <- renderPlot({
    grid::grid.newpage()
    grid::grid.draw(plot_abundance_class_legend())
  })
  

  
  ### OVERLAPS
  
  # Caching overlaps to speed up the loading of the tabs
  recall_base <- reactive({
    
    df <- switch(
      as.character(input$threshold_overlap_id),
      "60" = data_list$recall_fnr60,
      "70" = data_list$recall_fnr70,
      "80" = data_list$recall_fnr80,
      data_list$recall_fnr
    )
    
    df %>%
      filter(
        new_level %in% input$overlap_genes,
        tool_ref %in% input$tool_overlap,
        tool_comp %in% input$tool_overlap_calc
      ) %>%
      mutate(
        new_level = factor(new_level, levels = input$overlap_genes),
        tool_ref  = factor(tool_ref, levels = tools_levels),
        texture   = ifelse(tool_comp %in% tools_texture, "yes", "no")
      )
      
  }) %>%
    bindCache(
      input$threshold_overlap_id,
      input$overlap_genes,
      input$tool_overlap,
      input$tool_overlap_calc
    )
  
  
  plot_overlap_cstc  <- reactive({
    recall_base() %>%
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
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$overlap_cstc <- renderPlot({plot_overlap_cstc ()})
  output$overlap_cstc_overview <- renderPlot({plot_overlap_cstc ()})
  
  
  
  plot_overlap_cstc_summary  <- reactive({
    
    recall_base() %>%
      filter(!is.na(recall)) %>% 
      ungroup() %>% 
      group_by(tool_ref, new_level) %>% 
      summarise(recall = median(recall)) %>%
      ggplot(aes(x = "Class medians", fill = tool_ref, y = recall)) +
      geom_violin() +
      #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
      geom_jitter(size = 1.5, width = 0.4, height = 0) + 
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
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$overlap_cstc_summary <- renderPlot({plot_overlap_cstc_summary ()})
  output$overlap_cstc_summary_overview <- renderPlot({plot_overlap_cstc_summary ()})
  
  
  
  plot_overlap_csno  <- reactive({
    
    recall_base() %>%
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
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$overlap_csno <- renderPlot({plot_overlap_csno ()})
  output$overlap_csno_overview <- renderPlot({plot_overlap_csno ()})
  
  
  
  plot_overlap_csno_summary  <- reactive({
    
    recall_base() %>%
      filter(!is.na(recall)) %>% 
      ungroup() %>% 
      group_by(tool_ref, new_level) %>% 
      summarise(fnr = median(fnr)) %>%
      ggplot(aes(x = "Class medians", fill = tool_ref, y = fnr)) +
      geom_violin() +
      geom_jitter(size = 1.5, width = 0.4, height = 0) + 
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap], 
                        labels = tools_levels[tools_levels %in% input$tool_overlap]) +
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
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$overlap_csno_summary <- renderPlot({plot_overlap_csno_summary ()})
  output$overlap_csno_summary_overview <- renderPlot({plot_overlap_csno_summary ()})
  

  
  plot_overlap_legend <- reactive({
    
    tools_order <- match(input$tool_overlap_calc, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    
    grid.draw(g_legend(recall_base() %>%
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
  
  ## Adding this so it can be rendered separately well in the submenus
  output$overlap_legend <- renderPlot({plot_overlap_legend ()})
  
  plot_cstc_legend <- reactive({
    
    tools_order <- match(input$tool_overlap_calc, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    
    grid.draw(g_legend(recall_base() %>%
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
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$plot_cstc_legend <- renderPlot({plot_cstc_legend ()})
  
  
  ## Legend for CSNO
  plot_csno_legend <- reactive({
    
    tools_order <- match(input$tool_overlap_calc, tools_levels)
    shape_tools <- rep(21, length(tools_levels))
    shape_tools[tools_levels %in% tools_texture] <- 24
    shape_tools <- shape_tools[tools_order]
    
    grid.draw(g_legend(recall_base() %>%
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
  
  ## Adding this so it can be rendered seperately well in the submenus
  output$plot_csno_legend <- renderPlot({plot_csno_legend ()})
  
  
  ### Pan-/Core-resistome proportion
  pan_reactive <- reactive({
    
    if (input$threshold_pan_core_id == 60.0) {
      bind_rows(
        data_list$pan %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND")),
        data_list$pan60)
    } else if (input$threshold_pan_core_id == 70.0) {
      bind_rows(
        data_list$pan %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND")),
        data_list$pan70)
    } else if (input$threshold_pan_core_id == 80.0) {
      bind_rows(
        data_list$pan %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND")),
        data_list$pan80)
    } else {
      data_list$pan
    }
  })
  
  sumcore_reactive <- reactive({
    
    core_data <-
      if (input$threshold_pan_core_id == 60.0) {
        sum_core_adjust(bind_rows(
          data_list$core %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND")),
          data_list$core60
        ), input$threshold_samples, input$threshold_proportion )
      } else if (input$threshold_pan_core_id == 70.0) {
        sum_core_adjust(bind_rows(
          data_list$core %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND")),
          data_list$core70), input$threshold_samples, input$threshold_proportion
        )
      } else if (input$threshold_pan_core_id == 80.0) {
        sum_core_adjust(bind_rows(
          data_list$core %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND")),
          data_list$core80), input$threshold_samples, input$threshold_proportion)
      } else {
        sum_core_adjust(data_list$core, 
                        input$threshold_samples, 
                        input$threshold_proportion)
      }
  })
  
  
  output$plot_pan_core_proportion <- renderPlot({
    
    
    alluvial_pan_core_env(
      sumcore = sumcore_reactive() %>%
        filter(tool %in% input$tool_pan_core,
               habitat %in% input$single_environment_pan_core) %>%
        mutate(
          tool = factor(as.character(tool), 
                        levels = tools_levels[tools_levels %in% input$tool_pan_core])),
      pan     = pan_reactive() %>%
        filter(tool %in% input$tool_pan_core,
               habitat %in% input$single_environment_pan_core) %>%
        mutate(
          tool = factor(as.character(tool), 
                        levels = tools_levels[tools_levels %in% input$tool_pan_core])),
      h       = input$single_environment_pan_core,
      tools   = input$tool_pan_core,
      pal     = pal_10_complete,
      general_size = general_size,
      data_list$levels_unigenes) 
    
    
  })
  
  
}