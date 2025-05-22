# Load required libraries
library(shiny)
library(shinydashboard)
library(stringr)
library(Biostrings)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(dplyr)
library(viridis)

# Source helper functions
source("R/dna_functions.R")
source("R/visualization_helpers.R")

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "🧬 GenomeScope"),    
  
  dashboardSidebar(                             
    sidebarMenu(
      menuItem("Input", tabName = "input", icon = icon("dna")),
      menuItem("Basic Analysis", tabName = "basic", icon = icon("chart-bar")),
      menuItem("Advanced Analysis", tabName = "advanced", icon = icon("microscope")),
      menuItem("Structural Analysis", tabName = "structure", icon = icon("cube")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Input Tab
      tabItem(tabName = "input",
              box(width = 12,
                  textAreaInput("sequence", "Enter DNA Sequence:", 
                                rows = 5, 
                                placeholder = "Enter DNA sequence (A, T, G, C only)"),
                  actionButton("analyze_btn", "Analyze", 
                               icon = icon("play"),
                               class = "btn-primary"),
                  actionButton("example_btn", "Load Example", 
                               icon = icon("file-import"))
              ),
              box(width = 12,
                  verbatimTextOutput("validation_message"))
      ),
      
      # Basic Analysis Tab
      tabItem(tabName = "basic",
              fluidRow(
                box(width = 6,
                    plotlyOutput("nucleotide_plot"),
                    title = "Nucleotide Frequencies"
                ),
                box(width = 6,
                    plotOutput("gc_pie_plot"),
                    title = "GC Content Distribution"
                )
              ),
              fluidRow(
                box(width = 6,
                    DT::dataTableOutput("basic_stats"),
                    title = "Basic Statistics"
                ),
                box(width = 6,
                    DT::dataTableOutput("codon_usage_table"),
                    title = "Codon Usage"
                )
              )
      ),
      
      # Advanced Analysis Tab
      tabItem(tabName = "advanced",
              fluidRow(
                box(width = 12,
                    plotlyOutput("property_distribution"),
                    title = "Nucleotide Distribution"
                )
              ),
              fluidRow(
                box(width = 6,
                    plotlyOutput("gc_skew"),
                    title = "GC Skew"
                ),
                box(width = 6,
                    plotlyOutput("complexity_plot"),
                    title = "Sequence Complexity"
                )
              ),
              fluidRow(
                box(width = 12,
                    title = "Repeating Patterns",
                    DT::dataTableOutput("patterns_table"))
              )
      ),
      
      # Structure Tab
      tabItem(tabName = "structure",
              fluidRow(
                box(width = 12, height = "600px",
                    plotlyOutput("structure_3d"),
                    title = "3D DNA Structure"
                )
              ),
              fluidRow(
                box(width = 6,
                    plotOutput("circular_plot"),
                    title = "Circular DNA Visualization"
                ),
                box(width = 6,
                    plotOutput("codon_plot"),
                    title = "Codon Usage"
                )
              ),
              fluidRow(
                box(width = 12,
                    plotlyOutput("property_plot"),
                    title = "DNA Properties Overview"
                )
              ),
              fluidRow(
                box(width = 6,
                    title = "DNA Statistics",
                    DT::dataTableOutput("structure_stats")
                ),
                box(width = 6,
                    title = "Repeating Patterns",
                    DT::dataTableOutput("structure_patterns")
                )
              )
      ),
      
      # About Tab
      tabItem(tabName = "about",
              box(width = 12,
                  includeMarkdown("www/about.md"))
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Reactive values
  rv <- reactiveValues(
    sequence = NULL,
    analysis_results = NULL,
    advanced_results = NULL
  )
  
  # Load example sequence
  observeEvent(input$example_btn, {
    example_seq <- "ATGCATGCATGCATGCGAATTCGGATCCAAGCTTGCGGCCGCATGCAT"
    updateTextAreaInput(session, "sequence", value = example_seq)
  })
  
  # Analyze button click
  observeEvent(input$analyze_btn, {
    seq <- toupper(str_trim(input$sequence))
    
    # Validate sequence
    if (nchar(seq) == 0) {
      output$validation_message <- renderText("Please enter a DNA sequence")
      return()
    }
    
    if (!is_valid_dna(seq)) {
      output$validation_message <- renderText(
        "Invalid DNA sequence. Please enter only A, T, G, and C."
      )
      return()
    }
    
    # Store sequence and results
    rv$sequence <- seq
    rv$analysis_results <- analyze_dna(seq)
    rv$advanced_results <- list(
      properties = calculate_dna_properties(seq),
      structure = predict_dna_structure(seq),
      patterns = find_repeating_patterns(seq)
    )
    
    output$validation_message <- renderText("Analysis complete!")
  })
  
  # Basic Analysis Outputs
  output$nucleotide_plot <- renderPlotly({
    req(rv$analysis_results)
    counts <- rv$analysis_results$nucleotide_counts
    df <- data.frame(
      Nucleotide = names(counts),
      Count = as.numeric(counts)
    )
    
    plot_ly(df, x = ~Nucleotide, y = ~Count, type = "bar", 
            marker = list(color = c("#FF9999", "#FFFF99", "#9999FF", "#99FF99"))) %>%
      layout(title = "Nucleotide Frequencies",
             xaxis = list(title = "Nucleotide"),
             yaxis = list(title = "Count"))
  })
  
  output$gc_pie_plot <- renderPlot({
    req(rv$analysis_results)
    gc_content <- rv$analysis_results$gc_content
    at_content <- 100 - gc_content
    
    data <- data.frame(
      content = c("GC", "AT"),
      percentage = c(gc_content, at_content)
    )
    
    ggplot(data, aes(x = "", y = percentage, fill = content)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = c("GC" = "#9999FF", "AT" = "#FF9999")) +
      theme_void() +
      labs(title = "GC/AT Content") +
      theme(legend.position = "bottom")
  })
  
  output$basic_stats <- DT::renderDataTable({
    req(rv$analysis_results)
    data.frame(
      Property = c("Sequence Length", 
                   "GC Content", 
                   "Melting Temperature",
                   "Sequence Complexity"),
      Value = c(nchar(rv$sequence),
                paste0(round(rv$analysis_results$gc_content, 2), "%"),
                paste0(round(rv$advanced_results$properties$melting_temp, 2), "°C"),
                round(rv$advanced_results$properties$complexity, 4))
    )
  })
  
  output$codon_usage_table <- DT::renderDataTable({
    req(rv$sequence)
    codons <- substring(rv$sequence,
                        seq(1, nchar(rv$sequence)-2, by=3),
                        seq(3, nchar(rv$sequence), by=3))
    codon_counts <- table(codons)
    codon_freq <- prop.table(codon_counts) * 100
    
    data.frame(
      Codon = names(codon_counts),
      Count = as.vector(codon_counts),
      Frequency = sprintf("%.2f%%", as.vector(codon_freq))
    ) %>%
      arrange(desc(Count)) %>%
      DT::datatable(options = list(
        pageLength = 5,
        scrollX = TRUE
      ))
  })
  
  # Advanced Analysis Outputs
  output$property_distribution <- renderPlotly({
    req(rv$sequence)
    plots <- create_advanced_plots(rv$sequence)
    ggplotly(plots$distribution)
  })
  
  output$gc_skew <- renderPlotly({
    req(rv$sequence)
    plots <- create_advanced_plots(rv$sequence)
    ggplotly(plots$gc_skew)
  })
  
  output$complexity_plot <- renderPlotly({
    req(rv$sequence)
    # Calculate sequence complexity in sliding windows
    window_size <- 10
    n <- nchar(rv$sequence) - window_size + 1
    complexity_scores <- numeric(n)
    
    for(i in 1:n) {
      subseq <- substr(rv$sequence, i, i + window_size - 1)
      complexity_scores[i] <- calculate_sequence_complexity(subseq)
    }
    
    plot_data <- data.frame(
      Position = 1:n,
      Complexity = complexity_scores
    )
    
    plot_ly(plot_data, x = ~Position, y = ~Complexity,
            type = "scatter", mode = "lines",
            line = list(color = '#2ca02c')) %>%
      layout(title = "Sequence Complexity",
             xaxis = list(title = "Position"),
             yaxis = list(title = "Complexity Score"))
  })
  
  output$patterns_table <- DT::renderDataTable({
    req(rv$advanced_results)
    rv$advanced_results$patterns
  })
  
  # Structure Outputs
  output$structure_3d <- renderPlotly({
    req(rv$sequence)
    create_3d_structure(rv$sequence)
  })
  
  output$circular_plot <- renderPlot({
    req(rv$sequence)
    create_circular_plot(rv$sequence)
  })
  
  output$codon_plot <- renderPlot({
    req(rv$sequence)
    # Get codon usage data
    codons <- substring(rv$sequence,
                        seq(1, nchar(rv$sequence)-2, by=3),
                        seq(3, nchar(rv$sequence), by=3))
    codon_counts <- table(codons)
    codon_freq <- prop.table(codon_counts) * 100
    
    data.frame(
      Codon = factor(names(codon_counts), 
                     levels = names(codon_counts)[order(codon_freq, decreasing = TRUE)]),
      Frequency = as.vector(codon_freq)
    ) %>%
      ggplot(aes(x = Codon, y = Frequency, fill = Frequency)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis_c() +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Codon Usage Distribution",
           x = "Codon",
           y = "Frequency (%)")
  })
  
  output$property_plot <- renderPlotly({
    req(rv$sequence)
    
    # Window analysis
    window_size <- 10
    n <- nchar(rv$sequence) - window_size + 1
    
    properties <- data.frame(
      position = 1:n,
      gc = numeric(n),
      flex = numeric(n),
      tm = numeric(n)
    )
    
    for(i in 1:n) {
      subseq <- substr(rv$sequence, i, i + window_size - 1)
      properties$gc[i] <- calculate_gc_content(subseq)
      properties$tm[i] <- calculate_melting_temp(subseq)
      properties$flex[i] <- mean(predict_dna_flexibility(subseq)$flexibility)
    }
    
    plot_ly() %>%
      add_trace(data = properties,
                x = ~position,
                y = ~gc,
                name = "GC Content",
                type = "scatter",
                mode = "lines",
                line = list(color = '#1f77b4')) %>%
      add_trace(data = properties,
                x = ~position,
                y = ~tm/10,  # Scale temperature for visualization
                name = "Melting Temp",
                type = "scatter",
                mode = "lines",
                line = list(color = '#ff7f0e')) %>%
      add_trace(data = properties,
                x = ~position,
                y = ~flex * 100,  # Scale flexibility for visualization
                name = "Flexibility",
                type = "scatter",
                mode = "lines",
                line = list(color = '#2ca02c')) %>%
      layout(title = "DNA Properties Along Sequence",
             xaxis = list(title = "Position"),
             yaxis = list(title = "Relative Values"),
             showlegend = TRUE)
  })
  
  output$structure_stats <- DT::renderDataTable({
    req(rv$sequence)
    stats <- data.frame(
      Property = c(
        "Total Length",
        "GC Content",
        "AT Content",
        "Molecular Weight",
        "Melting Temperature",
        "Stability Score"
      ),
      Value = c(
        paste(nchar(rv$sequence), "bp"),
        paste0(round(rv$analysis_results$gc_content, 2), "%"),
        paste0(round(100 - rv$analysis_results$gc_content, 2), "%"),
        paste0(round(rv$advanced_results$properties$molecular_weight, 2), " Da"),
        paste0(round(rv$advanced_results$properties$melting_temp, 2), "°C"),
        round(rv$advanced_results$structure$stability, 3)
      )
    )
    DT::datatable(stats, options = list(dom = 't', ordering = FALSE))
  })
  
  output$structure_patterns <- DT::renderDataTable({
    req(rv$sequence)
    patterns <- find_repeating_patterns(rv$sequence) %>%
      arrange(desc(count)) %>%
      head(10)  # Show top 10 patterns
    
    DT::datatable(patterns,
                  options = list(
                    pageLength = 5,
                    scrollX = TRUE
                  )
    )
  })
}

# Run the app
shinyApp(ui = ui, server = server)
