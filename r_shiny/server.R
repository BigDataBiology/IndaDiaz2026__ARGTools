library(shiny)

# Detect the working directory
setwd("r_shiny") 

# Set the paths to data
DATA_DIR <- "code_R_analysis/output_abundance_diversity_resistome"
METADATA_PATH <- "data/metadata_GMGC10.sample.meta.tsv"

# Defining constants
general_size <- 10

# Label formatting function
lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}

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


# Defining the Server logic for the argCompare app

function(input, output, session) {
  #Loading data reactively
  # Main data loading reactive
  data <- reactive({
    # Show loading message
    showModal(modalDialog(
      title = "Loading Data...",
      "Please wait while the data is being loaded. This may take a moment.",
      footer = NULL,
      easyClose = FALSE
    ))
    
    # Validate paths exist
    if (!dir.exists(DATA_DIR)) {
      removeModal()
      showModal(modalDialog(
        title = "Error: Data Directory Not Found",
        paste0("Cannot find data directory: ", DATA_DIR, 
               "\n\nPlease check the BASE_DIR and DATA_DIR settings in ui.R"),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return(NULL)
    }
    
    # Try to load data with error handling
    tryCatch({
      # Load main datasets
      lst <- readRDS(file.path(DATA_DIR, "results_tools.rds"))
      metadata <- read.delim(METADATA_PATH)
      
      # Load the aggregated abundance per ARO and gene class
      abundance <- bind_rows(
        readRDS(file.path(DATA_DIR, "abundance_diversity_part1.rds")),
        readRDS(file.path(DATA_DIR, "abundance_diversity_part2.rds")),
        readRDS(file.path(DATA_DIR, "abundance_diversity_part3.rds"))
      )
      
      # Load core and pan resistome data
      core <- readRDS(file.path(DATA_DIR, "core_resistome.rds"))
      pan <- readRDS(file.path(DATA_DIR, "pan_resistome.rds"))
      
      # Process abundance data
      abundance <- abundance %>% 
        mutate(unigenes = distinct_unigenes_rarefied) %>%
        mutate(habitat = factor(habitat, levels = EN),
               habitat2 = factor(SO[habitat], levels = h2),
               tool = factor(tool, levels = tools_levels)) %>%
        mutate(location = ifelse(habitat2 %in% c("humans","mammals","wastewater","built-environment"), 
                                 "human-related","external")) %>%
        mutate(location = factor(location, levels = c("human-related","external"))) %>%
        filter(tool %in% tools_levels & !habitat %in% not_env) %>%
        filter(aggregation %in% "new_level") %>%
        mutate(habitat2 = factor(as.character(habitat2), h2))
      
      # Process core data
      core <- core %>% 
        rename(new_level = new_level_centroid, X = centroid) %>%
        filter(tool %in% tools_levels, !habitat %in% not_env) %>%
        mutate(habitat = factor(habitat, levels = EN2),
               tool = factor(tool, levels = tools_levels))
      
      # Process pan data
      pan <- pan %>% 
        filter(tool %in% tools_levels, !habitat %in% not_env,
               aggregation %in% "new_level_centroid") %>%
        mutate(habitat = factor(habitat, levels = EN2),
               tool = factor(tool, levels = tools_levels))
      
      removeModal()
      
      list(
        lst = lst,
        metadata = metadata,
        abundance = abundance,
        core = core,
        pan = pan
      )
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error Loading Data",
        paste0("Error: ", e$message, 
               "\n\nPlease check:\n",
               "1. All required .rds files are in: ", DATA_DIR, "\n",
               "2. Metadata file exists at: ", METADATA_PATH, "\n",
               "3. helper.R is in the same directory as ui.R and server.R"),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return(NULL)
    })
  })
  

  # FILTERED DATA REACTIVES
  
  # Filtered abundance data based on user selections
  filtered_abundance <- reactive({
    req(data())
    data()$abundance %>%
      filter(tool %in% input$select_tools,
             habitat %in% input$select_habitats)
  })
  
  # Filtered core resistome data
  filtered_core <- reactive({
    req(data())
    data()$core %>%
      filter(tool %in% input$select_tools,
             habitat %in% input$select_habitats)
  })
  
  # Filtered pan resistome data
  filtered_pan <- reactive({
    req(data())
    data()$pan %>%
      filter(tool %in% input$select_tools,
             habitat %in% input$select_habitats)
  })
  
  
  # Rendering the overview tab as valueboxes
  
  output$total_samples_box <- renderValueBox({
    req(data())
    n_samples <- n_distinct(data()$metadata$sample_id)
    valueBox(
      value = formatC(n_samples, format = "d", big.mark = ","),
      subtitle = "Total Samples",
      icon = icon("vial"),
      color = "blue"
    )
  })
  
  output$total_tools_box <- renderValueBox({
    req(input$select_tools)
    valueBox(
      value = length(input$select_tools),
      subtitle = "Selected Tools",
      icon = icon("wrench"),
      color = "green"
    )
  })
  
  output$total_habitats_box <- renderValueBox({
    req(input$select_habitats)
    valueBox(
      value = length(input$select_habitats),
      subtitle = "Selected Habitats",
      icon = icon("earth-americas"),
      color = "yellow"
    )
  })
  
  output$total_genes_box <- renderValueBox({
    req(filtered_abundance())
    n_genes <- filtered_abundance() %>%
      filter(unigenes > 0) %>%
      pull(unigenes) %>%
      sum()
    valueBox(
      value = formatC(n_genes, format = "d", big.mark = ","),
      subtitle = "Total Unique Genes",
      icon = icon("dna"),
      color = "red"
    )
  })
  
  # Rendering the Abundance Plot
  
  # Abundance plot
  output$plot_abundance <- renderPlot({
    req(filtered_abundance())
    
    abundance_tool_sample <- filtered_abundance() %>%
      group_by(tool, sample, habitat, habitat2) %>%
      summarise(normed10m = sum(normed10m), unigenes = sum(unigenes), .groups = "drop") %>%
      ungroup()
    
    ggplot(abundance_tool_sample, aes(x = habitat, y = normed10m, fill = tool)) +
      stat_summary(fun.data = calc_boxplot_stat, 
                   geom = "boxplot", 
                   position = position_dodge2(preserve = "single")) +
      scale_fill_manual(values = pal_10_q, labels = lab_fn(tools_levels)) +
      scale_y_continuous(labels = scales::label_number()) +
      labs(x = "", y = "Relative Abundance (per 10M reads)", fill = "") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = general_size),
        axis.text.y = element_text(size = general_size),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        legend.text = element_text(size = general_size),
        panel.border = element_rect(colour = "black", fill = NA)
      )
  })
  
  # Rendering the Diversity Plot
  output$plot_diversity <- renderPlot({
    req(filtered_abundance())
    
    abundance_tool_sample <- filtered_abundance() %>%
      group_by(tool, sample, habitat, habitat2) %>%
      summarise(normed10m = sum(normed10m), unigenes = sum(unigenes), .groups = "drop") %>%
      ungroup()
    
    ggplot(abundance_tool_sample, aes(x = habitat, y = unigenes, fill = tool)) +
      stat_summary(fun.data = calc_boxplot_stat, 
                   geom = "boxplot", 
                   position = position_dodge2(preserve = "single")) +
      scale_fill_manual(values = pal_10_q, labels = lab_fn(tools_levels)) +
      labs(x = "", y = "Alpha Diversity (Number of Unique Genes)", fill = "") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = general_size),
        axis.text.y = element_text(size = general_size),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        legend.text = element_text(size = general_size),
        panel.border = element_rect(colour = "black", fill = NA)
      )
  })
  

  # The Pan & Core Resistome Plot
  
  # Pan vs Core plot
  output$plot_pancore <- renderPlot({
    req(filtered_pan(), filtered_core())
    
    # Calculate core summary
    sumcore <- sum_core_adjust(filtered_core(), input$min_samples, input$core_threshold)
    
    # Calculate pan summary
    sumpan2 <- filtered_pan() %>%
      ungroup() %>%
      group_by(tool, habitat, aggregation, epoch) %>%
      summarise(s = sum(unigenes), .groups = "drop") %>%
      ungroup() %>%
      group_by(tool, habitat, aggregation) %>%
      summarise(md = median(s), mn = mean(s), sd = sd(s), .groups = "drop")
    
    # Combine pan and core
    pan_core <- sumpan2 %>%
      left_join(sumcore %>%
                  ungroup() %>%
                  group_by(tool, habitat) %>%
                  summarise(core = sum(unigenes), .groups = "drop"),
                by = c("tool", "habitat")) %>%
      mutate(core = ifelse(is.na(core), 0, core)) %>%
      mutate(prop = core / md) %>%
      ungroup() %>%
      group_by(tool, habitat) %>%
      mutate(texture = ifelse(tool %in% tools_texture, "yes", "no"))
    
    ggplot(pan_core, aes(x = habitat)) +
      geom_col(aes(y = md, fill = tool), position = position_dodge(), alpha = 1) +
      geom_col(aes(y = core, fill = tool), position = position_dodge(), alpha = 0.4) +
      scale_fill_manual(values = pal_10_q) +
      labs(x = "", y = "Number of Genes", fill = "",
           title = "Pan-resistome (full bars) vs Core-resistome (faded bars)") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = general_size),
        axis.text.y = element_text(size = general_size),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        plot.title = element_text(size = general_size + 2, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA)
      )
  })
  
  # Core prevalence plot
  output$plot_core_prevalence <- renderPlot({
    req(filtered_core())
    
    core_data <- filtered_core() %>%
      filter(cut == input$core_threshold, cnt > input$min_samples) %>%
      group_by(tool, habitat) %>%
      summarise(n_genes = n_distinct(X), .groups = "drop")
    
    ggplot(core_data, aes(x = habitat, y = n_genes, fill = tool)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = pal_10_q) +
      labs(x = "", y = "Core Genes", fill = "") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = general_size),
        axis.text.y = element_text(size = general_size)
      )
  })
  
  # Pan growth plot
  output$plot_pan_growth <- renderPlot({
    req(filtered_pan())
    
    pan_summary <- filtered_pan() %>%
      group_by(tool, habitat, epoch) %>%
      summarise(genes = sum(unigenes), .groups = "drop") %>%
      group_by(tool, habitat) %>%
      summarise(mean_genes = mean(genes), .groups = "drop")
    
    ggplot(pan_summary, aes(x = habitat, y = mean_genes, fill = tool)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = pal_10_q) +
      labs(x = "", y = "Mean Pan-resistome Size", fill = "") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = general_size),
        axis.text.y = element_text(size = general_size)
      )
  })
  
  # TOOL OVERLAP PLOTS
  
  # Overlap heatmap
  output$plot_overlap_heatmap <- renderPlot({
    req(data(), input$overlap_habitat)
    
    # Filter for selected habitat
    habitat_data <- data()$lst[[input$overlap_habitat]]
    
    if (!is.null(habitat_data)) {
      overlap_matrix <- return_overlap_tools(habitat_data)
      
      ggplot(overlap_matrix, aes(x = tool_ref, y = tool_comp, fill = recall)) +
        geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.2f", recall)), color = "white", size = 3) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5) +
        labs(x = "Reference Tool", y = "Comparison Tool", fill = "Recall") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = general_size),
          axis.text.y = element_text(size = general_size),
          axis.title = element_text(size = general_size + 1, face = "bold")
        )
    } else {
      plot.new()
      text(0.5, 0.5, "No data available for selected habitat", cex = 1.5)
    }
  })
  
  # Venn diagram
  output$plot_venn <- renderPlot({
    req(data(), input$overlap_habitat, input$venn_tools)
    
    if (length(input$venn_tools) < 2 || length(input$venn_tools) > 5) {
      plot.new()
      text(0.5, 0.5, "Please select 2-5 tools for Venn diagram", cex = 1.2)
      return()
    }
    
    habitat_data <- data()$lst[[input$overlap_habitat]]
    
    if (!is.null(habitat_data)) {
      gene_lists <- habitat_data %>%
        filter(tool %in% input$venn_tools) %>%
        group_by(tool) %>%
        summarise(genes = list(unique(query)), .groups = "drop") %>%
        deframe()
      
      if (length(gene_lists) >= 2) {
        ggVennDiagram(gene_lists, label_alpha = 0) +
          scale_fill_gradient(low = "white", high = "red") +
          theme(legend.position = "right")
      } else {
        plot.new()
        text(0.5, 0.5, "Not enough tools selected", cex = 1.2)
      }
    } else {
      plot.new()
      text(0.5, 0.5, "No data available for selected habitat", cex = 1.2)
    }
  })
  
  # Overlap statistics table
  output$table_overlap_stats <- renderReactable({
    req(data(), input$overlap_habitat)
    
    habitat_data <- data()$lst[[input$overlap_habitat]]
    
    if (!is.null(habitat_data)) {
      overlap_stats <- return_overlap_tools(habitat_data) %>%
        filter(tool_ref != tool_comp) %>%
        select(tool_ref, tool_comp, jaccard, recall) %>%
        mutate(across(where(is.numeric), ~round(., 3)))
      
      reactable(
        overlap_stats,
        columns = list(
          tool_ref = colDef(name = "Reference Tool", minWidth = 120),
          tool_comp = colDef(name = "Comparison Tool", minWidth = 120),
          jaccard = colDef(name = "Jaccard Index", minWidth = 100),
          recall = colDef(name = "Recall", minWidth = 80)
        ),
        searchable = TRUE,
        filterable = TRUE,
        defaultPageSize = 10,
        striped = TRUE,
        highlight = TRUE
      )
    }
  })
  

  # GENE CLASS PLOTS
  
  # Gene class plot
  output$plot_geneclass <- renderPlot({
    req(filtered_abundance(), input$geneclass_habitat)
    
    gene_summary <- filtered_abundance() %>%
      filter(habitat == input$geneclass_habitat) %>%
      group_by(tool, gene) %>%
      summarise(total_abundance = sum(normed10m),
                total_diversity = sum(unigenes), .groups = "drop") %>%
      arrange(desc(total_abundance)) %>%
      slice_head(n = 20)
    
    ggplot(gene_summary, aes(x = reorder(gene, total_abundance), 
                             y = total_abundance, fill = tool)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = pal_10_q) +
      coord_flip() +
      labs(x = "Gene Class", y = "Total Abundance", fill = "") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text = element_text(size = general_size),
        axis.title = element_text(size = general_size + 1, face = "bold")
      )
  })
  
  # Gene class table
  output$table_geneclass <- renderReactable({
    req(filtered_abundance(), input$geneclass_habitat)
    
    gene_table <- filtered_abundance() %>%
      filter(habitat == input$geneclass_habitat) %>%
      group_by(gene, tool) %>%
      summarise(
        samples = n_distinct(sample),
        mean_abundance = mean(normed10m),
        mean_diversity = mean(unigenes),
        .groups = "drop"
      ) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      arrange(desc(mean_abundance))
    
    reactable(
      gene_table,
      columns = list(
        gene = colDef(name = "Gene Class", minWidth = 200),
        tool = colDef(name = "Tool", minWidth = 150),
        samples = colDef(name = "Samples", minWidth = 80),
        mean_abundance = colDef(name = "Mean Abundance", minWidth = 120),
        mean_diversity = colDef(name = "Mean Diversity", minWidth = 120)
      ),
      searchable = TRUE,
      filterable = TRUE,
      defaultPageSize = 15,
      striped = TRUE,
      highlight = TRUE
    )
  })
  
  # SUMMARY TABLES
  
  # Sample summary table
  output$table_sample_summary <- renderReactable({
    req(data())
    
    sample_summary <- data()$metadata %>%
      filter(habitat %in% input$select_habitats) %>%
      group_by(habitat) %>%
      summarise(n_samples = n(), .groups = "drop") %>%
      arrange(desc(n_samples))
    
    reactable(
      sample_summary,
      columns = list(
        habitat = colDef(name = "Habitat", minWidth = 150),
        n_samples = colDef(name = "Number of Samples", minWidth = 120)
      ),
      striped = TRUE,
      highlight = TRUE
    )
  })
  
  # Tool summary table
  output$table_tool_summary <- renderReactable({
    req(filtered_abundance())
    
    tool_summary <- filtered_abundance() %>%
      group_by(tool) %>%
      summarise(
        n_samples = n_distinct(sample),
        n_genes = n_distinct(gene),
        total_abundance = sum(normed10m),
        mean_diversity = mean(unigenes),
        .groups = "drop"
      ) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      arrange(desc(total_abundance))
    
    reactable(
      tool_summary,
      columns = list(
        tool = colDef(name = "Tool", minWidth = 180),
        n_samples = colDef(name = "Samples", minWidth = 80),
        n_genes = colDef(name = "Unique Genes", minWidth = 100),
        total_abundance = colDef(name = "Total Abundance", minWidth = 120),
        mean_diversity = colDef(name = "Mean Diversity", minWidth = 120)
      ),
      striped = TRUE,
      highlight = TRUE
    )
  })
}

