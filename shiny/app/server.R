

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
  
  #ARGs TAB
  output$plot_count_genes_tool <- renderPlot({
    unigenes %>%
      filter(tool %in% input$tools_unigenes) %>% 
      ungroup() %>% 
      group_by(tools_labels, texture, tools_db, tool2) %>% 
      summarise(n = n_distinct(query)) %>%
      ggplot(aes(x = tool2, y = n, fill = tools_db, pattern = texture)) +
      geom_col_pattern(position = position_dodge2(preserve = "single", width = 0.8), 
                       width = 0.8, pattern_color = "black", pattern_fill = pattern_fill, 
                       pattern_size = 0.12, color = "black") +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      facet_grid(. ~ tools_db, scales = "free_x", space = "free") +
      scale_fill_manual(values = pal_7) +
      xlab("Pipelines") + ylab("Number of ARGs") + 
      scale_y_continuous(expand = c(0.01, 0.01), 
                         breaks = c(25000,50000,75000,100000,125000), 
                         labels = scales::comma) + 
      guides(fill = guide_legend(
    override.aes = list(
      pattern = rep("none", 7),
      fill  = pal_7)), pattern = "none") + 
  theme1 +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5))
  }) %>% bindCache(input$tools_unigenes)
  
  
  output$plot_gene_class_proportion <- renderPlot({
    
    req(input$tools_unigenes)
    req(input$gene_classes_filter)
    plot_data <- unigenes_propotion %>%
      filter(tool %in% input$tools_unigenes) %>%
      filter(new_level %in% input$gene_classes_filter)

    req(nrow(plot_data) > 0)
   
    plot_data %>%
      ggplot(aes(x = tool2, y = new_level, fill = p)) +
      geom_tile(color = "grey") +
      scale_fill_gradientn(
        colors = brewer.pal(9, "YlOrBr"),
        labels = percent_format(accuracy = 1),
        breaks = c(0.0001, 0.2, 0.4, 0.6)
      ) +
      facet_grid(. ~ tools_db, scales = "free_x", space = "free") +
      scale_y_discrete(labels = function(x) {
        out <- ifelse(
          x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
          paste0("italic('", x, "')"),
          ifelse(x %in% "abcF", "ABC-F", paste0("'", x, "'"))
        )
        parse(text = out)
      }) +
      geom_text(
        data = plot_data %>% filter(p >= 0.03),
        aes(label = scales::percent(p, accuracy = 1), color = p >= 0.40),
        size = general_size / .pt,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c("black", "white")) +
      labs(fill = "") +
      ylab("Proportion of ARG class") +
      xlab("Pipelines") +
      theme(
        text = element_text(size = general_size, color = "black"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size, face = "bold"),
        strip.text = element_text(size = general_size, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = general_size, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 9, hjust = 1),
        plot.margin = margin(10, 20, 10, 10, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(2, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0)
      )
    
  }, height = 1000, res = 96) %>% bindCache(input$tools_unigenes, input$gene_classes_filter)
  
  
  #ABUNDANCE TAB
  filtered_abundance_data <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    abundance_tool_sample %>%
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance)
  })
  
  filtered_class_abundance_data <- reactive({
    req(input$tool_abundance, input$environment_abundance, input$abundance_genes)
    
    abundance_class_summary %>% 
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance) %>% 
      filter(gene %in% input$abundance_genes) %>% 
      ungroup() %>%
      mutate(tools_labels = droplevels(tools_labels))
  })
  
  output$plot_abundance <- renderPlot({
    
    ggplot(filtered_abundance_data(), aes(x = tool2, y = abundance, fill = tools_db, pattern = texture)) + 
      geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 1, padding = 0), 
                           pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                           pattern_spacing = 0.2,
                           pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA,
                           linewidth = 0.15) +
      scale_x_discrete(expand = expansion(add = 1)) +
      
      geom_jitter(data = filtered_abundance_data() %>% 
                    ungroup() %>%
                    group_by(habitat, tool) %>% 
                    filter(abundance < quantile(abundance, 0.75) + 1.5*IQR(abundance)) %>%
                    group_modify(~ dplyr::slice_sample(.x, n = min(100, nrow(.x)))), 
                  width = 0.35, size = 0.4, alpha = 0.4) +
      facet_grid(habitat ~ tools_db  , scales = "free") +
      scale_y_continuous(labels = scales::comma) + 
      scale_fill_manual(values = pal_7) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      xlab("Pipeline") +
      ylab("Relative abundance (per million)") +
      theme1 + 
    theme(panel.border = element_blank(),
          strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt"))
  }, height = 600, res = 96) %>% bindCache(input$tool_abundance, input$environment_abundance)
  
  
  output$plot_abundance_gene_class <- renderPlot({
    
    ggplot(filtered_class_abundance_data(), aes(y = fct_rev(tool2), fill = tools_db, pattern = texture)) +
      geom_boxplot_pattern( 
        aes(xmin = q25, xlower = q25, xmiddle = q50, xupper = q75, xmax = q75, pattern = texture),
        stat = "identity", 
        position = position_dodge2(preserve = "single", padding = 0), 
        pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.07,
        pattern_density = 0.15,
        pattern_size =  0.07, color = "black", outliers = FALSE, outlier.shape = NA,
        linewidth = 0.1) +
      scale_fill_manual(values = pal_7) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      facet_grid(gene ~ habitat, scales = "free", space = "free",
                 labeller = labeller(
                   gene = as_labeller(function(x) {
                     out <- ifelse(
                       x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim", "vat", "mph", "qnr"),
                       paste0("italic('", x, "')"),
                       ifelse(x == "abcF", "ABC-F",
                              paste0("'", x, "'"))
                     )
                     return(out)
                   }, default = label_parsed)
                 )) +
      xlab("Relative abundance (per million)") +
      ylab("") + 
      labs(fill = "") + 
      scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
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
    
  }, height = 600, res = 96) %>% bindCache(input$tool_abundance, input$environment_abundance, input$abundance_genes)
  
  #PAN- & CORE- RESISTOME TAB
  output$pan_core <- renderPlot({
    
    pan_core <- sumpan2 %>% 
      filter(tool %in% input$tool_pan_core,
             habitat %in% input$environment_pan_core) %>% 
      left_join((sumcore2 %>% filter(tool %in% input$tool_pan_core,
                                     cnt == 450, 
                                     cut == 0.5,
                                     habitat %in% input$environment_pan_core) %>%
                   ungroup() %>%
                   select(-c(cut, cnt, tools_labels, tools_db, texture))), by = c("tool", "habitat", "tool2")) %>%
      mutate(core = ifelse(is.na(core), 0, core)) %>% 
      mutate(prop = core / md) %>%
      ungroup() 
    
    pan_core_df_plot <- pan_core %>% 
      select(!c(md, sd)) %>% 
      pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
      mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
      mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
      mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
      filter(value > 0)
    
    
    p4a1 <- pan_core_df_plot %>% filter(metric %in% "Pan-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x =  value)) +
      geom_point_rast(aes(fill = tool2, shape = texture),  color = "black", stroke = 0.3, size = 3) + 
      facet_grid(habitat ~ metric, scales = "free") +
      scale_shape_manual(values = c("no" = 21, "yes" = 24, 'y70' = 22, 'y80' = 23,  'y90' = 25)) +
      scale_fill_manual(values = pal_10_q_2[names(pal_10_q_2) %in% pan_core_df_plot$tool2]) +
      theme_minimal() +
      ylab("") +
      xlab("Number of ARGs") +
      labs(fill = "") +
      ggtitle("a") + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      theme_minimal() +
      guides(fill = guide_legend(
        override.aes = list(
          shape = shape_tools_2[names(shape_tools_2) %in% pan_core_df_plot$tool2],
          fill  = pal_10_q_2[names(pal_10_q_2) %in% pan_core_df_plot$tool2],
          color  = pal_10_q_2[names(pal_10_q_2) %in% pan_core_df_plot$tool2])), shape = "none") + 
      theme(
        legend.position = "none",
        text = element_text(size = 16, color = "black"),
        title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16 , face = "bold"),
        strip.text = element_text(size = 16, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        panel.border = element_rect(color =  "black"),
        plot.margin = margin(5, 5, 5, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank()) +
      theme(axis.text.y = element_blank(), 
            strip.text.x = element_text(size = 16, angle = 0, vjust = 0, hjust = 0.5),
            strip.text.y = element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0.5),
            legend.position = "bottom",
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.border = element_blank()) +
      scale_x_reverse(labels = label_comma()) 
    
    p4a2 <- pan_core_df_plot %>% filter(!metric %in% "Pan-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x =  value)) +
      geom_point_rast(aes(fill = tool2, shape = texture),  color = "black", stroke = 0.3, size = 3) + 
      facet_grid(habitat ~ metric, scales = "free") +
      scale_fill_manual(values = pal_10_q_2[names(pal_10_q_2) %in% pan_core_df_plot$tool2]) +
      scale_shape_manual(values = c("no" = 21, "yes" = 24, 'y70' = 22, 'y80' = 23,  'y90' = 25)) +
      theme_minimal() +
      ylab("") +
      xlab("Number of ARGs") +
      labs(fill = "") +
      ggtitle("b") + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      theme_minimal() +
      guides(fill = guide_legend(
        override.aes = list(
          shape = shape_tools_2[names(shape_tools_2) %in% pan_core_df_plot$tool2],
          fill  = pal_10_q_2[names(pal_10_q_2) %in% pan_core_df_plot$tool2],
          color  = pal_10_q_2[names(pal_10_q_2) %in% pan_core_df_plot$tool2])), shape = "none") + 
      theme(
        legend.position = "none",
        text = element_text(size = 16, color = "black"),
        title = element_text(size = 16 + 2, face = "bold"),
        axis.title = element_text(size = 16 , face = "bold"),
        strip.text = element_text(size = 16, angle = 0),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 16, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 16),
        panel.border = element_rect(color =  "black"),
        plot.margin = margin(5, 5, 5, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.minor.y = element_blank()) +
      theme(axis.text.y = element_blank(), 
            strip.text.x = element_text(size = 16, angle = 0, vjust = 0, hjust = 0.5),
            strip.text.y = element_blank(),
            legend.position = "none",
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.border = element_blank())
    
    
    panel1 <- (p4a1 + theme(legend.position = "none") | p4a2) + patchwork::plot_layout(widths = c(1,1))
    p4 <- (panel1 / (g_legend(p4a1) )) + patchwork::plot_layout(heights = c(15,1))
    p4
    
  }) %>% bindCache(input$tool_pan_core, input$environment_pan_core)
  
  #OVERLAP TAB
  filtered_overlap_data <- reactive({
    req(input$tool_overlap, input$overlap_genes) 
    
    recall_fnr %>% 
      filter(tool_ref %in% input$tool_overlap, tool_comp %in% basic_tools) %>%
      filter(new_level %in% input$overlap_genes)
  })
  
  output$overlap_gene_class <- renderPlot({
    
    ggplot(filtered_overlap_data(), 
           aes(x = recall*100, y = new_level)) + 
      geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture),
                           position = position_dodge2(preserve = "single"),
                           color = "black", outliers = T,
                           pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.07,
                           pattern_density = 0.15,
                           pattern_size =  0.07,
                           linewidth = 0.1) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 
                                      'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) + 
      facet_grid(tool_ref2 ~ " ", scales = "free_y") +
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap]) +
      scale_y_discrete(drop = FALSE) +
      xlab("Percentage (%)") +
      ylab("ARG class") + 
      theme_minimal() +
      theme5 +
      theme( panel.grid = element_blank(),
             strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
             strip.text.y = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 0))
  }, height = 600, res = 96) %>% bindCache(input$tool_overlap, input$overlap_genes)
  
  
  output$overlap <- renderPlot({
    
    ggplot(filtered_overlap_data(), # 
           aes(x = recall*100, y = fct_rev(tool_comp2))) +
      geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture),
                           position = position_dodge2(preserve = "single"),
                           color = "black", outliers = T,
                           pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.07,
                           pattern_density = 0.15,
                           pattern_size =  0.07,
                           linewidth = 0.1) +
      facet_grid(tools_db_comp ~ tools_labels_ref, scales = "free_y", space = "free") +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 
                                      'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) + 
      scale_fill_manual(values = pal_10_q[tools_levels %in% input$tool_overlap]) +
      xlab("%") +
      ylab("Pipeline covered") +
      theme_minimal() +
      theme5 +
      theme(panel.grid = element_blank(),
            plot.margin = margin(0, 10, 0, 0, unit = "pt"),
            strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
            strip.text.y = element_text(size = general_size, angle=0, vjust = 0.5, hjust = 0))
    
  }, height = 600, res = 96) %>% bindCache(input$tool_overlap, input$overlap_genes)
  
}


