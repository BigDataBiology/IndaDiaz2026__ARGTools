

server <- function(input, output, session) {
  
  # Wait 500 milliseconds after the user stops clicking before updating plots
  genes_debounced <- reactive({ input$abundance_genes }) %>% debounce(500)
  tools_debounced <- reactive({ input$tool_abundance }) %>% debounce(500)
  
  # Helper to translate the UI threshold (e.g., "0.0" or "60.0") to the list name ("default" or "60")
  get_thr_name <- function(thr_val) {
    thr_str <- as.character(thr_val)
    thr_str <- gsub("\\.0$", "", thr_str) # This strips the trailing .0
    
    if (thr_str %in% c("0", "default")) {
      return("default")
    }
    return(thr_str)
  }
  
  ##
  
  output$plot_count_genes_tool <- renderPlot({
    p1 <- unigenes %>%
    filter(tool %in% input$tools_unigenes) %>% 
    ungroup() %>% 
    group_by(tools_labels, texture, tools_db) %>% 
    summarise(n = n_distinct(query)) %>%
    ggplot(aes (x = tools_labels, y = n, fill = tools_db, pattern = texture )) +
    geom_col_pattern(position = position_dodge2(preserve = "single", width = 0.8), 
                     width = 0.8, pattern_color = "black", pattern_fill = pattern_fill, 
                     pattern_size =  0.12, color = "black") +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
    facet_grid( . ~ tools_db, scales = "free_x", space = "free") +
    scale_fill_manual(values = pal_7) +
    xlab("") + ylab("Number of ARGs") + ggtitle("a") + 
    scale_y_continuous(labels = scales::comma) + 
    theme(
      legend.position = "none",
      text = element_text(size = general_size, color = "black"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size , face = "bold"),
      strip.text = element_text(size = general_size, angle = 0),
      panel.background = element_blank(),
      axis.text.x = element_text(size = general_size, angle = 90),
      axis.text.y = element_text(size = general_size),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.spacing = unit(0, "pt"),
      legend.text = element_text(size = general_size),
      panel.grid.minor.y = element_blank(), panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background.x = element_blank(),
          strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0)) 
  
    p2 <- unigenes_propotion %>% 
      filter(tool %in% input$tools_unigenes) %>% 
      ggplot(aes(x = tools_labels, y = new_level, fill = p)) + 
      geom_tile(color = "grey") + 
      scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                            labels = percent_format(accuracy = 1),
                            breaks = c(0.0001, 0.2, 0.4, 0.6)) + 
      facet_grid(.~tools_db, scales = "free_x", space = "free") +
      geom_text(
        data = unigenes_propotion %>% filter(tool %in% input$tools_unigenes) %>% filter(p >= 0.03),  # only these get labels
        aes(label = scales::percent(p, accuracy = 1),
            color = p >= 0.40), 
        size = general_size / .pt, show.legend = F) +
      scale_color_manual(values = c( "black", "white")) +
      labs(fill = "") + 
      ggtitle("b") + 
      ylab("Proportion of ARG class") + 
      xlab("") + 
      theme(
        text = element_text(size = general_size, color = "black"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size , face = "bold"),
        strip.text = element_text(size = general_size, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = general_size, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom", panel.grid = element_blank(),
        panel.border =  element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0))
    
    (p1 | p2) + patchwork::plot_layout(widths = c(1,1.5))
  })
  


  output$plot_abundance <- renderPlot({

    p1 <- abundance_tool_sample %>%
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance) %>% 
      ggplot(aes(x = tools_labels, y = abundance, fill = tools_db, pattern = texture)) + 
      geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 1, padding = 0), 
                           width = 1, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                           pattern_spacing = 0.2,
                           pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA,
                           linewidth = 0.15) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_jitter(data = abundance_tool_sample %>% filter(tool %in% input$tool_abundance) %>% 
                    filter(habitat %in% input$environment_abundance) %>% 
                    ungroup() %>%
                    group_by(habitat, tool) %>% 
                    filter(abundance < quantile(abundance, 0.75) + 1.5*IQR(abundance)) %>%
                    mutate(n_group = n()) %>%
                    group_modify(~ dplyr::slice_sample(.x, n = min(100, nrow(.x)))), 
                  width = 0.35, size = 0.4, alpha = 0.1) + 
      facet_grid(habitat ~ tools_db  , scales = "free") +
      scale_y_continuous(labels = scales::comma) + 
      scale_fill_manual(values = pal_7) +
      ggtitle("a") + 
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      theme_minimal() +
      theme(
        legend.position = "none",
        text = element_text(size = general_size, color = "black"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size , face = "bold"),
        strip.text = element_text(size = general_size, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = general_size, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size),
        panel.border = element_rect(color =  "black"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(10, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank()) +
      theme(panel.border = element_blank(),
            strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10, unit = "pt"))  
    

    
    p2 <- abundance_class_summary %>% 
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance) %>% 
      filter(gene %in% input$abundance_genes) %>% 
      ungroup() %>%
      mutate(tools_labels = droplevels(tools_labels)) %>% 
      ggplot(aes(y = fct_rev(tool), fill = tools_db, pattern = texture)) +
      geom_boxplot_pattern( 
        aes(xmin = q25, xlower = q25, xmiddle = q50, xupper = q75, xmax = q75, pattern = texture),
        stat = "identity", 
        position = position_dodge2(preserve = "single", width = 0.5, padding = 0), 
        width = 1, pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.07,
        pattern_density = 0.15,
        pattern_size =  0.07, color = "black", outliers = FALSE, outlier.shape = NA,
        linewidth = 0.1) +
      scale_fill_manual(values = pal_7) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      facet_grid(gene ~ habitat, scales = "free", space = "free") +
      xlab("Relative abundance") +
      ylab("") + 
      labs(fill = "") + 
      scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
      ggtitle("b") +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      theme_minimal() +
      theme(
        legend.position = "none",
        text = element_text(size = general_size, color = "black"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size , face = "bold"),
        strip.text = element_text(size = general_size, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = general_size, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size),
        panel.border = element_rect(color =  "black"),
        plot.margin = margin(10, 10, 10, 10, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(10, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank()) +
      theme(
        legend.position = "none",
        strip.text.y = element_text(angle = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border =  element_blank(),
        axis.text.y = element_blank())   
      
    (p1 | p2) + patchwork::plot_layout(widths = c(1,1)) 
    
    
  })
  
  
    
  output$pan_core <- renderPlot({
    
    sumcore <- sum_core_adjust(core %>% 
                                 filter(tool %in% input$tool_pan_core,
                                        habitat %in% input$environment_pan_core), 
                               input$threshold_samples, input$threshold_proportion) %>%
      mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
             texture = ifelse(tool %in% tools_texture, "yes", "no"),
             tools_db = factor(tools_db[tool], levels = tools_db_factor))
    
    pan_core <- sumpan2 %>% 
      filter(tool %in% input$tool_pan_core,
             habitat %in% input$environment_pan_core) %>% 
      left_join((sumcore %>% 
                   ungroup() %>% 
                   group_by(tool, habitat) %>% 
                   summarise(core = sum(unigenes))), by = c("tool", "habitat")) %>%
      mutate(core = ifelse(is.na(core), 0, core)) %>% 
      mutate(prop = core / md) %>%
      ungroup() 
    
    pan_core_df_plot <- pan_core %>% 
      select(!c(md, sd)) %>% 
      pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
      mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
      mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
      mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
      filter(value > 0) %>% 
      mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                              ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                     ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                            texture))))
    
    
    p4a1 <- pan_core_df_plot %>% filter(metric %in% "Pan-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x =  value)) +
      geom_point(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 3) + 
      facet_grid(habitat ~ metric, scales = "free") +
      scale_fill_manual(values = pal_10_q[names(pal_10_q) %in% pan_core_df_plot$tool],
                        labels = tools_labels[tools_levels %in% pan_core_df_plot$tool]) +
      scale_shape_manual(values = c("no" = 21, "yes" = 24, 'y70' = 22, 'y80' = 23,  'y90' = 25)) +
      theme_minimal() +
      ylab("") +
      xlab("Number of ARGs") +
      labs(fill = "") +
      ggtitle("a") + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      theme_minimal() +
      guides(fill = guide_legend(
        override.aes = list(
          shape = shape_tools[names(shape_tools) %in% pan_core_df_plot$tool],
          fill  = pal_10_q[names(pal_10_q) %in% pan_core_df_plot$tool],
          color  = pal_10_q[names(pal_10_q) %in% pan_core_df_plot$tool])), shape = "none") + 
      theme(
        legend.position = "none",
        text = element_text(size = general_size, color = "black"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size , face = "bold"),
        strip.text = element_text(size = general_size, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = general_size, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size),
        panel.border = element_rect(color =  "black"),
        plot.margin = margin(5, 5, 5, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank()) +
      theme(axis.text.y = element_blank(), 
            strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
            legend.position = "bottom",
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.border = element_blank()) +
      scale_x_reverse(labels = label_comma()) 
    
    p4a2 <- pan_core_df_plot %>% filter(!metric %in% "Pan-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x =  value)) +
      geom_point(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 3) + 
      facet_grid(habitat ~ metric, scales = "free") +
      scale_fill_manual(values = pal_10_q[names(pal_10_q) %in% pan_core_df_plot$tool],
                        labels = tools_labels[tools_levels %in% pan_core_df_plot$tool]) +
      scale_shape_manual(values = c("no" = 21, "yes" = 24, 'y70' = 22, 'y80' = 23,  'y90' = 25)) +
      theme_minimal() +
      ylab("") +
      xlab("Number of ARGs") +
      labs(fill = "") +
      ggtitle("a") + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      theme_minimal() +
      guides(fill = guide_legend(
        override.aes = list(
          shape = shape_tools[names(shape_tools) %in% pan_core_df_plot$tool],
          fill  = pal_10_q[names(pal_10_q) %in% pan_core_df_plot$tool],
          color  = pal_10_q[names(pal_10_q) %in% pan_core_df_plot$tool])), shape = "none") + 
      theme(
        legend.position = "none",
        text = element_text(size = general_size, color = "black"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size , face = "bold"),
        strip.text = element_text(size = general_size, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = general_size, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size),
        panel.border = element_rect(color =  "black"),
        plot.margin = margin(5, 5, 5, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank()) +
      theme(axis.text.y = element_blank(), 
            strip.text.x = element_blank(),
            legend.position = "none",
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.border = element_blank())
    
    
    panel1 <- (p4a1 + theme(legend.position = "none") | p4a2) + patchwork::plot_layout(widths = c(1,1))
    p4 <- (panel1 / (g_legend(p4a1) )) + patchwork::plot_layout(heights = c(15,1))
    p4
    
  })
  
    
}