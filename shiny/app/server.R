

server <- function(input, output, session) {
  
  
  #ARGs TAB
  count_genes_reactive <- reactive({
    req(input$tools_unigenes)
    unigenes %>%
      filter(tool %in% input$tools_unigenes) %>%
      ungroup() %>%
      group_by(tools_labels, texture, tools_db, tool2) %>%
      summarise(n = n_distinct(query), .groups = "drop")
  })
  
  output$plot_count_genes_tool <- renderPlot({
    count_genes_reactive() %>%
      ggplot(aes(x = tools_labels, y = n, fill = tools_db, pattern = texture)) +
      geom_col_pattern(position = position_dodge2(preserve = "single", width = 0.8), 
                       width = 0.8, pattern_color = "black", pattern_fill = pattern_fill, 
                       pattern_size = 0.12, color = "black") +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      facet_grid(. ~ tools_db, scales = "free_x", space = "free") +
      scale_fill_manual(values = pal_7, drop = TRUE) +
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
  
  output$download_count_genes <- downloadHandler(
    filename = function() paste0("ARG_counts_", Sys.Date(), ".csv"),
    content = function(file) {
      count_genes_reactive() %>%
        select(tools_labels, tools_db, n) %>%
        rename(
          "Tool"    = tools_labels,
          "Tools DB" = tools_db,
          "Number of ARGs"   = n
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  output$plot_gene_class_proportion <- renderPlot({
    
    req(input$tools_unigenes)
    req(input$gene_classes_filter)
    plot_data <- unigenes_proportion %>%
      filter(tool %in% input$tools_unigenes) %>%
      filter(new_level %in% input$gene_classes_filter)
    
    req(nrow(plot_data) > 0)
    
    plot_data %>%
      ggplot(aes(x = tools_labels, y = new_level, fill = p)) +
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
  
  output$download_gene_class_proportion <- downloadHandler(
    filename = function() paste0("gene_class_proportion_", Sys.Date(), ".csv"),
    content = function(file) {
      unigenes_proportion %>%
        filter(tool %in% input$tools_unigenes) %>%
        filter(new_level %in% input$gene_classes_filter) %>%
        ungroup() %>%
        select(tool, new_level, n, p) %>%
        rename(
          "Tool"                = tool,
          "ARG Class"            = new_level,
          "Number of ARGs"       = n,
          "ARG Class Proportion" = p
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  #ABUNDANCE TAB
  filtered_abundance_data <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    abundance_tool_sample %>%
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance) %>%
      group_by(tool, habitat) %>%
      mutate(N_samples = n_distinct(sample)) %>%
      ungroup() %>%
      mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(N_samples), ")"))

  })
  
  filtered_class_abundance_data <- reactive({
    req(input$tool_abundance, input$environment_abundance, input$abundance_genes)
    
    abundance_class_summary %>% 
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance) %>% 
      filter(gene %in% input$abundance_genes) %>% 
      ungroup() %>%
      mutate(tools_labels = droplevels(tools_labels)) %>%
      left_join(
        abundance_tool_sample %>%
          group_by(habitat) %>%
          summarise(N_samples = n_distinct(sample)),
        by = "habitat"
      ) %>%
      mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(N_samples), ")"))
  })
  
  output$plot_abundance <- renderPlot({
    
    ggplot(filtered_abundance_data(), aes(x = tools_labels, y = abundance, fill = tools_db, pattern = texture)) + 
      geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 1, padding = 0), 
                           pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                           pattern_spacing = 0.2,
                           pattern_size =  0.3, color = "black", outliers = FALSE, 
                           linewidth = 0.15) +
      scale_x_discrete(expand = expansion(add = 1)) +
      
      geom_jitter(data = filtered_abundance_data() %>% 
                    ungroup() %>%
                    group_by(habitat, tool) %>% 
                    filter(abundance < quantile(abundance, 0.75) + 1.5*IQR(abundance)) %>%
                    group_modify(~ dplyr::slice_sample(.x, n = min(100, nrow(.x)))), 
                  width = 0.35, size = 0.4, alpha = 0.4) +
      facet_grid(habitat_label ~ tools_db, scales = "free") +
      scale_y_continuous(labels = scales::comma) + 
      scale_fill_manual(values = pal_7) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      xlab("Pipelines") +
      ylab("Relative abundance\n(aligned reads per million)") +
      theme1 + 
      theme(panel.border = element_blank(),
            strip.text.x = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 0),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt"))
    
  }, height = 600, res = 96) %>% bindCache(input$tool_abundance, input$environment_abundance)
  
  output$download_abundance <- downloadHandler(
    filename = function() paste0("abundance_", Sys.Date(), ".csv"),
    content = function(file) {
      filtered_abundance_data() %>%
        select (sample, tool, abundance, richness, habitat) %>%
        rename(
          "Sample ID"    = sample,
          "Tool" = tool,
          "Abundance (aligned reads per million)" = abundance,
          "Richness" = richness,
          "Environment" = habitat
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  output$plot_abundance_gene_class <- renderPlot({
    
    ggplot(filtered_class_abundance_data(), aes(y = fct_rev(tools_labels), fill = tools_db, pattern = texture)) +
      geom_boxplot_pattern( 
        aes(xmin = q25, xlower = q25, xmiddle = q50, xupper = q75, xmax = q75, pattern = texture),
        stat = "identity", 
        position = position_dodge2(preserve = "single", padding = 0), 
        pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.07,
        pattern_density = 0.15,
        pattern_size =  0.07, color = "black", outliers = FALSE,
        linewidth = 0.1) +
      scale_fill_manual(values = pal_7) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      facet_grid(gene ~ habitat_label, scales = "free", space = "free",
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
      xlab("Relative abundance\n(aligned reads per million)") +
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
  
  output$download_class_abundance <- downloadHandler(
    filename = function() paste0("gene_class_abundance_", Sys.Date(), ".csv"),
    content = function(file) {
      filtered_class_abundance_data() %>%
        ungroup() %>%
        select(tool, habitat, gene, q50, q25, q75, w1, w2) %>%
        rename(
          "Tool"            = tool,
          "Environment"      = habitat,
          "ARG"              = gene,
          "Median (q50)"     = q50,
          "Q25"              = q25,
          "Q75"              = q75,
          "Whisker Low (w1)" = w1,
          "Whisker High (w2)"= w2
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  #PAN- & CORE- RESISTOME TAB
  pan_core_reactive <- reactive({
    req(input$tool_pan_core, input$environment_pan_core, input$threshold_samples)
    
    pan_core <- sumpan2 %>% 
      filter(tool %in% input$tool_pan_core,
             habitat %in% input$environment_pan_core) %>% 
      left_join((sumcore2 %>% filter(tool %in% input$tool_pan_core,
                                     cnt %in% input$threshold_samples, 
                                     cut == 0.5,
                                     habitat %in% input$environment_pan_core) %>%
                   ungroup() %>%
                   select(-c(cut, cnt, tools_labels, tools_db, texture))), by = c("tool", "habitat", "tool2")) %>%
      mutate(core = ifelse(is.na(core), 0, core),
             prop = core / md) %>%
      ungroup() 
    
    pan_core %>% 
      select(!c(md, sd)) %>% 
      pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
      mutate(metric = case_when(
        metric == "mn"   ~ "Pan-resistome",
        metric == "core" ~ "Core-resistome"
      )) %>%
      mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
      filter(value > 0)
  })
  
  
  output$pan_core <- renderPlot({
  
    
    pan_core_df_plot <- pan_core_reactive()
      
    req(nrow(pan_core_df_plot) > 0)
    
    present_tools <- tools_levels[tools_levels %in% unique(pan_core_df_plot$tool)]
    legend_fills  <- pal_10_q[present_tools]
    legend_shapes <- shape_tools[present_tools]
    legend_labels <- as.vector(tools_labels[present_tools])
    
    theme_pan_core <- theme_minimal() +
      theme(
        text          = element_text(size = 16, color = "black"),
        title         = element_text(size = 16, face = "bold"),
        axis.title    = element_text(size = 16, face = "bold"),
        axis.text.x   = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y   = element_blank(),
        strip.text.x  = element_text(size = 16, angle = 0, vjust = 0, hjust = 0.5),
        panel.background = element_blank(),
        panel.border     = element_blank(),
        panel.spacing    = unit(0, "pt"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin      = margin(5, 5, 5, 5, unit = "pt"),
        legend.position  = "none",
        legend.text      = element_text(size = general_size),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin    = margin(0, 0, 0, 0, unit = "pt")
      )
    
    scale_fill_pan  <- scale_fill_manual(values = pal_10_q_2[levels(droplevels(pan_core_df_plot$tool2))])
    scale_shape_pan <- scale_shape_manual(values = c("no" = 21, "yes" = 24, "y70" = 22, "y80" = 23, "y90" = 25))
    
    p4a1 <- pan_core_df_plot %>% filter(metric == "Pan-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x = value)) +
      geom_point_rast(aes(fill = tool2, shape = texture), color = "black", stroke = 0.3, size = 3) + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      facet_grid(habitat ~ metric, scales = "free") +
      scale_fill_pan + scale_shape_pan +
      scale_x_reverse(labels = label_comma()) +
      labs(y = "", x = "Number of ARGs", fill = "", title = "a") +
      theme_pan_core +
      theme(strip.text.y = element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0.5))
    
    p4a2 <- pan_core_df_plot %>% filter(metric == "Core-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x = value)) +
      geom_point_rast(aes(fill = tool2, shape = texture), color = "black", stroke = 0.3, size = 3) + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      facet_grid(habitat ~ metric, scales = "free") +
      scale_fill_pan + scale_shape_pan +
      scale_x_continuous(labels = label_comma()) +
      labs(y = "", x = "Number of ARGs", fill = "", title = "b",
     subtitle = paste0("Core-resistome threshold: >", input$threshold_samples, " subsamples")) +
      theme_pan_core +
      theme(strip.text.y = element_blank())
    
    legend_df <- data.frame(
      tool2 = factor(legend_labels, levels = legend_labels),
      x = 1, y = 1,
      texture = ifelse(legend_shapes == 21, "no",
                       ifelse(legend_shapes == 24, "yes",
                              ifelse(legend_shapes == 22, "y70",
                                     ifelse(legend_shapes == 23, "y80", "y90"))))
    )
    
    legend_dummy <- ggplot(legend_df, aes(x, y, fill = tool2, shape = texture)) +
      geom_point(size = 3, color = "black", stroke = 0.3) +
      scale_fill_manual(values = setNames(legend_fills, legend_labels)) +
      scale_shape_manual(values = c("no" = 21, "yes" = 24, "y70" = 22, "y80" = 23, "y90" = 25)) +
      guides(fill = guide_legend(
        override.aes = list(shape = legend_shapes, fill = legend_fills, color = legend_fills),
        label.theme = element_text(size = general_size),
        nrow = 2
      ), shape = "none") +
      labs(fill = "") +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.text      = element_text(size = general_size),
            legend.key.height = unit(0.35, "cm"),
            legend.key.width  = unit(0.35, "cm"))
    
    panel1 <- (p4a1 | p4a2) + patchwork::plot_layout(widths = c(1, 1))
    (panel1 / g_legend(legend_dummy)) + patchwork::plot_layout(heights = c(15, 1))
    
  }) 
  
  output$download_pan_resistome <- downloadHandler(
    filename = function() paste0("pan_resistome_", Sys.Date(), ".csv"),
    content = function(file) {
      pan_core_reactive() %>%
        filter(metric == "Pan-resistome") %>%
        select(tool, habitat, prop, metric, value) %>%
        rename(
          "Tool"        = tool,
          "Environment" = habitat,
          "Proportion"  = prop,
          "Metric"      = metric,
          "Mean richness across the subsamples"  = value
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  output$download_core_resistome <- downloadHandler(
    filename = function() paste0("core_resistome_", Sys.Date(), ".csv"),
    content = function(file) {
      pan_core_reactive() %>%
        filter(metric == "Core-resistome") %>%
        select(tool, habitat, prop, metric, value) %>%
        rename(
          "Tool"        = tool,
          "Environment" = habitat,
          "Proportion"  = prop,
          "Metric"      = metric,
          "Number of ARGs in the subsample" = value
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  #OVERLAP TAB
  filtered_overlap_data <- reactive({
    req(input$tool_overlap, input$overlap_genes) 
    
    csc_fnr %>% 
      filter(tool_ref %in% input$tool_overlap, tool_comp %in% basic_tools) %>%
      filter(new_level %in% input$overlap_genes)
  })
  
  output$overlap <- renderPlot({
    overlap_data <- filtered_overlap_data() %>%
      group_by(tool_ref, tools_labels_ref, tools_db_comp) %>%
      mutate(n_obs = n()) %>%
      mutate(n_obs = paste0('n = ', n_obs)) %>%
      ungroup()
    
    ggplot(overlap_data,  
           aes(x = csc*100, y = fct_rev(tools_labels_comp))) +
      geom_boxplot(aes(fill = tool_ref, pattern = texture),
                   position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                   width = 0.7, color = "black", outliers = FALSE) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe',
                                      'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      geom_text(aes(x = 50, y = 1, label = n_obs), size = 4) +
      facet_grid(tools_db_comp ~ tools_labels_ref, scales = "free_y", space = "free") +
      scale_fill_manual(values = pal_10_q) +
      scale_y_discrete(drop = T) +
      xlab("Percentage (%)") +
      ylab("Pipeline covered") +
      theme_minimal() +
      theme5 +
      theme(panel.grid = element_blank(),
            plot.margin = margin(0, 10, 0, 0, unit = "pt"),
            strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
            strip.text.y = element_text(size = general_size, angle=90, vjust = 0.5, hjust = 0.5))
    
  }, height = 600, res = 96) %>% bindCache(input$tool_overlap, input$overlap_genes)
  
  output$download_overlap <- downloadHandler(
    filename = function() paste0("overlap_csc_", Sys.Date(), ".csv"),
    content = function(file) {
      filtered_overlap_data() %>%
        select (new_level, tool_ref, tool_comp, csc, ref_n_class, comp_n_class, ref_n_all, comp_n_all) %>%
        rename (
          "ARG Class" = new_level,
          "Reference Tool" = tool_ref,
          "Compared Tool" = tool_comp,
          "Class-Specific Coverage" = csc,
          "Reference N (Class)" = ref_n_class,
          "Compared N (Class)" = comp_n_class,
          "Reference N (All)" = ref_n_all,
          "Compared N (All)" = comp_n_all
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )

  
  output$overlap_gene_class <- renderPlot({
    
    overlap_gene <- filtered_overlap_data() %>%
      group_by(tool_ref, tools_labels_ref, new_level) %>%
      mutate(n_obs = n()) %>%
      mutate(n_obs = paste0('n = ', n_obs)) %>%
      ungroup()
    
    
    ggplot(overlap_gene, aes(x = csc*100, y = new_level)) + 

      geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture),
                           position = position_dodge2(preserve = "single", width = 0.3, padding = 0),
                           width = 0.7, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                           pattern_spacing = 0.2,
                           pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +

      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe',
                                      'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      geom_text(aes(x = 50, y = 1, label = n_obs), size = general_size / .pt) + 
      facet_grid(new_level ~ tools_labels_ref, scales = "free_y", space = "free",
                 labeller = labeller(
                   new_level = as_labeller(function(x) {
                     out <- ifelse(
                       x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim", "vat", "mph", "qnr"),
                       paste0("italic('", x, "')"),
                       ifelse(x == "abcF", "ABC-F",
                              paste0("'", x, "'"))
                     )
                     return(out)
                   }, default = label_parsed)
                 )) +
      scale_fill_manual(values = pal_10_q)  +
      xlab("Percentage (%)") +
      ylab("ARG class") + 
      theme_minimal() +
      theme5 +
      theme( panel.grid = element_blank(),
             strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5 , angle = 0),
             strip.text.y = element_text(size = general_size, vjust = 0.5, hjust = 0, angle = 0),
             panel.spacing = unit(5, "pt"),
             axis.text.y = element_blank())
      
  }, height = 600, res = 96) %>% bindCache(input$tool_overlap, input$overlap_genes)
  
  output$download_overlap_gene_class <- downloadHandler(
    filename = function() paste0("overlap_gene_class_", Sys.Date(), ".csv"),
    content = function(file) {
      filtered_overlap_data() %>%
        select(new_level, tool_ref, csc) %>%
        rename (
          "ARG Class" = new_level,
          "Reference Tool" = tool_ref,
          "Class-Specific Coverage" = csc
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  
}


