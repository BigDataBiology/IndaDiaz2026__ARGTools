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
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe', 'y70' = 'circle', 'y80' = 'circle', 'y90' = 'circle')) +
      facet_grid(. ~ tools_db, scales = "free_x", space = "free") +
      scale_fill_manual(values = pal_7, drop = TRUE) +
      xlab("Pipelines") + ylab("Number of ARGs") + 
      scale_y_continuous(expand = c(0.01, 0.01), 
                         labels = scales::comma) + 
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
      filter(new_level %in% input$gene_classes_filter)%>%
      ungroup() %>%
      mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
      mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
      mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
      mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
      mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
      mutate(new_level = gsub("target-modifying enzyme","target-modif.\nenzyme", new_level)) %>%
      mutate(new_level = gsub("self-resistance","self-resistance", new_level)) 
    
    req(nrow(plot_data) > 0)
    
    plot_data %>%
      ggplot(aes(x = tools_labels, y = new_level, fill = p)) +
      geom_tile(color = "grey") +
      scale_fill_gradientn(
        colors = brewer.pal(9, "YlOrBr"),
        labels = percent_format(accuracy = 1),
        breaks = c(0.0001, 0.2, 0.4, 0.6)) +
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
      theme1 +
      theme(
        axis.text.y = element_text(size = general_size - 1),
        plot.margin = margin(0, 0, 0, 20, unit = "pt"),
        legend.position = "bottom",
        legend.key.height = unit(0.6, "cm"),
        legend.key.width  = unit(0.7, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5)
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
          "Tool" = tool,
          "ARG Class" = new_level,
          "Number of ARGs" = n,
          "ARG Class Proportion" = p
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  output$jaccard_plot <- renderPlot({

    JI_all_plot %>%
      filter(as.numeric(tool_lab_ref) < as.numeric(tool_lab_comp)) %>%
      ggplot(aes(x = tool_lab_ref, y = tool_lab_comp, fill = jaccard)) +
      geom_tile(color = "grey") +
      scale_fill_gradientn(
        colors = brewer.pal(9, "YlOrBr"),
        labels = percent_format(accuracy = 1),
        breaks = c(0, 0.2, 0.4, 0.6, 0.8),
      ) +
      geom_text(
        data = JI_all_plot %>%
          filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools,
                 as.numeric(tool_lab_ref) < as.numeric(tool_lab_comp),
                 jaccard >= 0.05),
        aes(label = scales::percent(jaccard, accuracy = 1),
            color  = jaccard >= 0.50),
        size = general_size / .pt, show.legend = FALSE
      ) +
      scale_color_manual(values = c("black", "white")) +
      scale_y_discrete(labels = function(x) {
        x <- gsub("ABRicate-\n", "ABRicate-", x)
        x <- gsub("-\nPlus", "Plus", x)
        x
      }, expand = 0) +
      labs(fill = "") +
      xlab("") + ylab("") +
      theme1 +
      theme(legend.position = "bottom",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width  = unit(0.7, "cm"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            plot.margin = margin(0, 0, 0, 10, unit = "pt"))
  })
  
  output$download_jaccard_plot <- downloadHandler(
    filename = function() paste0("jaccard_index_", Sys.Date(), ".csv"),
    content = function(file) {
      JI_all_plot %>%
        filter(as.numeric(tool_lab_ref) < as.numeric(tool_lab_comp)) %>%
        select(tool_ref, tool_comp, jaccard) %>%
        rename(
          "Pipeline A" = tool_ref,
          "Pipeline B" = tool_comp,
          "Jaccard Index" = jaccard
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
  }) %>% bindCache(input$tool_abundance, input$environment_abundance)
  
  filtered_abundance_jitter <- reactive({
    filtered_abundance_data() %>%
      ungroup() %>%
      group_by(habitat, tool) %>%
      filter(abundance < quantile(abundance, 0.75) + 1.5 * IQR(abundance)) %>%
      group_modify(~ slice_sample(.x, n = min(100, nrow(.x)))) %>%
      ungroup()
  }) %>% bindCache(input$tool_abundance, input$environment_abundance)
  
  filtered_class_abundance_data <- reactive({
    req(input$tool_abundance, input$environment_abundance, input$abundance_genes)
    
    dynamic_gene_levels <- input$abundance_genes
    dynamic_gene_levels <- gsub(" beta-lactamase","", dynamic_gene_levels)
    dynamic_gene_levels <- gsub("cell wall ","cell\nwall\n", dynamic_gene_levels)
    dynamic_gene_levels <- gsub("MFS efflux pump","MFS\nefflux", dynamic_gene_levels)
    dynamic_gene_levels <- gsub("efflux pump","efflux", dynamic_gene_levels)
                                
    abundance_class_summary %>% 
      filter(tool %in% input$tool_abundance) %>% 
      filter(habitat %in% input$environment_abundance) %>% 
      filter(gene %in% input$abundance_genes) %>% 
      mutate(gene = gsub(" beta-lactamase","", gene)) %>%
      mutate(gene = gsub("cell wall ","cell\nwall\n", gene)) %>%
      mutate(gene = gsub("MFS efflux pump","MFS\nefflux", gene)) %>%
      mutate(gene = gsub("efflux pump","efflux", gene)) %>%
      mutate(gene = factor(gene, levels = dynamic_gene_levels)) %>%
      ungroup() %>%
      mutate(tools_labels = droplevels(tools_labels)) %>%
      left_join(habitat_n_samples, by = "habitat") %>%
      mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(N_samples), ")")) 
    
  }) %>% bindCache(input$tool_abundance, input$environment_abundance, input$abundance_genes)
  
  output$plot_abundance <- renderPlot({
    
    ggplot(filtered_abundance_data(), aes(x = tools_labels, y = abundance, fill = tools_db, pattern = texture)) + 
      geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 1, padding = 0), 
                           pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                           pattern_spacing = 0.2,
                           pattern_size =  0.3, color = "black", outliers = FALSE, 
                           linewidth = 0.15) +
      scale_x_discrete(expand = expansion(add = 1)) +
      
      geom_jitter(data = filtered_abundance_jitter(),
                  width = 0.35, size = 0.4, alpha = 0.4) +
      facet_grid(habitat_label ~ tools_db, scales = "free") +
      scale_y_continuous(labels = scales::comma) + 
      scale_fill_manual(values = pal_7) +
      scale_pattern_shared +
      xlab("Pipelines") +
      ylab("Relative abundance\n(aligned reads per million)") +
      theme1 + 
      theme(panel.border = element_blank(),
            strip.text.x = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 0.5),
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
          "Sample ID" = sample,
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
                       x %in% c("rpoB", "van", "fos", "erm", "cat", 
                                "aph", "ant", "aac", "lnu", "nim", 
                                "vat", "mph", "qnr"),
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
        plot.margin = margin(10, 10, 10, 10, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(10, "pt"),
        legend.text = element_text(size = general_size),
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
          "Tool" = tool,
          "Environment" = habitat,
          "ARG" = gene,
          "Median (q50)" = q50,
          "Q25" = q25,
          "Q75" = q75,
          "Whisker Low (w1)" = w1,
          "Whisker High (w2)"= w2
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  #RICHNESS MENU
  filtered_richness_data <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    abundance_tool_sample %>%
      filter(tool %in% input$tool_abundance) %>%
      filter(habitat %in% input$environment_abundance) %>%
      group_by(habitat) %>%
      mutate(N_samples = n_distinct(sample)) %>%
      ungroup() %>%
      mutate(habitat = as.character(habitat)) %>%
      mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(N_samples), ")")) %>%
      group_by(tool, habitat, habitat_label, tools_db, tools_labels, texture) %>%
      summarise(
        median = ifelse(quantile(richness, 0.5)  < 0, 0, quantile(richness, 0.5)),
        q25 = ifelse(quantile(richness, 0.25) < 0, 0, quantile(richness, 0.25)),
        q75 = ifelse(quantile(richness, 0.75) < 0, 0, quantile(richness, 0.75)),
        w1 = ifelse(quantile(richness, 0.25) - 1.5*IQR(richness) < 0, 0,
                        quantile(richness, 0.25) - 1.5*IQR(richness)),
        w2 = ifelse(quantile(richness, 0.75) + 1.5*IQR(richness) < 0, 0,
                        quantile(richness, 0.75) + 1.5*IQR(richness)),
        .groups = "drop"
      )
    
  }) %>% bindCache(input$tool_abundance, input$environment_abundance)
  
  filtered_richness_jitter <- reactive({
    req(input$tool_abundance, input$environment_abundance)
    
    habitat_labels <- filtered_richness_data() %>%
      distinct(habitat, habitat_label)
    
    abundance_tool_sample %>%
      filter(tool %in% input$tool_abundance) %>%
      filter(habitat %in% input$environment_abundance) %>%
      mutate(habitat = as.character(habitat)) %>%
      group_by(habitat, tool) %>%
      filter(richness < quantile(richness, 0.75) + 1.5 * IQR(richness)) %>%
      ungroup() %>%
      group_by(habitat, tool) %>%
      slice_sample(n = 100) %>%
      ungroup() %>%
      left_join(habitat_labels, by = "habitat")
    
  }) %>% bindCache(input$tool_abundance, input$environment_abundance)
  
  output$plot_richness <- renderPlot({
    
    ggplot(filtered_richness_data(), aes(x = tools_labels, y = median, fill = tools_db)) +
      geom_boxplot_pattern(
        aes(ymin = w1, lower = q25, middle = median, upper = q75, ymax = w2, pattern = texture),
        stat = "identity",
        position = position_dodge2(preserve = "single", padding = 0),
        pattern_color = "black", pattern_fill = "black",
        pattern_density = 0.000000001, pattern_spacing = 0.2,
        pattern_size = 0.3, color = "black", outliers = FALSE,
        outlier.shape = NA, linewidth = 0.15
      ) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_jitter(
        data = filtered_richness_jitter(),
        aes(y = richness),
        width = 0.35, size = 0.4, alpha = 0.4
      ) +
      facet_grid(habitat_label ~ tools_db, scales = "free", space = "free_x") +
      scale_y_continuous(labels = scales::comma) +
      scale_fill_manual(values = pal_7) +
      scale_pattern_shared +
      xlab("Pipelines") + ylab("Richness") +
      theme1 +
      theme(
        panel.border = element_blank(),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
      )
    
  })
  
  output$download_richness <- downloadHandler(
    filename = function() paste0("richness_", Sys.Date(), ".csv"),
    content = function(file) {
      filtered_richness_data() %>%
        select(tool, habitat, median, q25, q75, w1, w2) %>%
        rename(
          "Pipeline" = tool,
          "Habitat" = habitat,
          "Median" = median,
          "Q25" = q25,
          "Q75" = q75,
          "Whisker Low" = w1,
          "Whisker High" = w2
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
    
  
  
  #PAN- & CORE- RESISTOME TAB
  pan_core_reactive <- reactive({
    req(input$tool_pan_core, input$environment_pan_core, input$threshold_samples, input$threshold_proportion)
    
    pan_core <- sumpan2 %>% 
      filter(tool %in% input$tool_pan_core,
             habitat %in% input$environment_pan_core) %>% 
      left_join((sumcore2 %>% filter(tool %in% input$tool_pan_core,
                                     cnt %in% input$threshold_samples, 
                                     cut %in% input$threshold_proportion,
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
    
  }) %>% bindCache(input$tool_pan_core, input$environment_pan_core, input$threshold_samples, input$threshold_proportion)
  
  
  output$pan_core <- renderPlot({
  
    
    pan_core_df_plot <- pan_core_reactive()
    
    req(nrow(pan_core_df_plot) > 0)
    
    present_tools <- tools_levels[tools_levels %in% unique(pan_core_df_plot$tool)]
    legend_fills  <- pal_10_q[present_tools]
    legend_shapes <- shape_tools[present_tools]
    legend_labels <- as.vector(tools_labels[present_tools])
    scale_fill_pan  <- scale_fill_manual(values = pal_10_q_2[levels(droplevels(pan_core_df_plot$tool2))])
    scale_shape_pan <- scale_shape_manual(values = c("no" = 21, "yes" = 24, "y70" = 22, "y80" = 23, "y90" = 25))
    
    p4a1 <- pan_core_df_plot %>% filter(metric == "Pan-resistome") %>% 
      ggplot(aes(y = fct_rev(tool), x = value)) +
      geom_point_rast(aes(fill = tool2, shape = texture), color = "black", stroke = 0.3, size = 3) + 
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      facet_grid(habitat ~ metric, scales = "free") +
      geom_segment(aes(x = 0, xend = value, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
      scale_color_manual(values = scale_fill_pan) +
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
      geom_segment(aes(x = 0, xend = value, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
      scale_color_manual(values = scale_fill_pan) +
      scale_fill_pan + scale_shape_pan +
      scale_x_continuous(labels = label_comma()) +
      labs(y = "", x = "Number of ARGs", fill = "", title = "b") +
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
            legend.text = element_text(size = general_size),
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
          "Tool"  = tool,
          "Environment" = habitat,
          "Proportion"  = prop,
          "Metric"  = metric,
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
          "Tool" = tool,
          "Environment" = habitat,
          "Proportion"  = prop,
          "Metric" = metric,
          "Number of ARGs in the subsample" = value
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )
  
  
  #OVERLAP TAB
  filtered_overlap_data <- reactive({
    req(input$tool_overlap, input$overlap_genes, input$tool_overlap_comp) 
    
    csc_fnr %>% 
      filter(tool_ref %in% input$tool_overlap, tool_comp %in% input$tool_overlap_comp) %>%
      filter(new_level %in% input$overlap_genes)%>%
      mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
      mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
      mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
      mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
      mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
      mutate(new_level = gsub("target-modifying enzyme","target-modif.\nenzyme", new_level)) %>%
      mutate(new_level = gsub("self-resistance","self-resistance", new_level))
    
  }) %>% bindCache(input$tool_overlap, input$overlap_genes, input$tool_overlap_comp)
  
  
  output$overlap_gene_class <- renderPlot({

    overlap_gene <- filtered_overlap_data() %>%
      group_by(tool_ref, new_level) %>%
      mutate(n_obs = n_distinct(tool_comp)) %>%
      mutate(n_obs = paste0('n = ', n_obs)) %>%
      ungroup()

    label_data <- overlap_gene %>%
      group_by(tool_ref, tools_labels_ref, new_level) %>%
      summarise(
        n_obs  = paste0("n = ", n()),
        x_pos  = quantile(csc * 100, 0.75, na.rm = TRUE) + 3,
        .groups = "drop"
      )

    ggplot(overlap_gene, aes(x = csc*100, y = new_level)) +

      geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture),
                           position = position_dodge2(preserve = "single", width = 0.3, padding = 0),
                           width = 0.7, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                           pattern_spacing = 0.2,
                           pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +

      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe',
                                      'y70' = 'crosshatch', 'y80' = 'crosshatch',  'y90' = 'crosshatch')) +
      geom_text(data = label_data, aes(x = x_pos, y = 1, label = n_obs), hjust = 0,
                size = 3.5, inherit.aes = FALSE) +
      scale_x_continuous(limits = c(0, 115), breaks = c(0, 25, 50, 75, 100)) +
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

  }, height = 600, res = 96) %>% bindCache(input$tool_overlap, input$overlap_genes, input$tool_overlap_comp)

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
  
  #SUPPLEMENTARY TABLES TAB
  output$table_s1 <- renderReactable({
    reactable(
      table_s1,
      searchable = TRUE,
      filterable = TRUE,
      striped    = TRUE,
      highlight  = TRUE,
      defaultPageSize = 20,
      columns = list(
        Sample  = colDef(name = "Sample"),
        Habitat = colDef(name = "Habitat")
      )
    )
  })
  
  output$download_table_s1 <- downloadHandler(
    filename = function() paste0("TableS1_samples_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(table_s1, file, row.names = FALSE)
  )
  
  output$table_s2 <- renderReactable({
    reactable(
      table_s2,
      searchable = TRUE,
      filterable = TRUE,
      striped    = TRUE,
      highlight  = TRUE,
      defaultPageSize = 20,
      columns = list(
        gene_class   = colDef(name = "Gene Class"),
        ARO_Term_ID  = colDef(name = "ARO Term ID"),
        abbreviation = colDef(name = "Abbreviation")
      )
    )
  })
  
  output$download_table_s2 <- downloadHandler(
    filename = function() paste0("TableS2_ARG_classes_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(table_s2, file, row.names = FALSE)
  )
  
  output$table_s3 <- renderReactable({
    
    table_s3_base %>%
      rename(
        "ARG (Unigenes)" = query,
        "Gene Class"    = new_level,
        "Pipeline"      = tool
      ) %>%
      
      reactable(
        searchable      = TRUE,
        filterable      = TRUE,
        striped         = TRUE,
        highlight       = TRUE,
        defaultPageSize = 20
      )
  })
  
  output$download_table_s3 <- downloadHandler(
    filename = function() paste0("ARGs_per_pipeline_", Sys.Date(), ".csv"),
    content  = function(file) {
      table_s3_base %>%
        rename(
          "ARG (Unigene)" = query,
          "Gene Class"    = new_level,
          "Pipeline"      = tool
        ) %>%
        write.csv(file, row.names = FALSE)
    }
  )

}


