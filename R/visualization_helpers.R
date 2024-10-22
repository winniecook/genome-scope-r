#' Visualization Helper Functions
#' Additional visualization functions for DNA analysis

library(ggplot2)
library(plotly)
library(viridis)
library(tidyr)

#' Create DNA property heatmap
#' @param sequence DNA sequence
#' @return ggplot object
create_property_heatmap <- function(sequence) {
  # Calculate properties in sliding windows
  window_size <- 10
  windows <- substring(sequence, 
                      seq_len(nchar(sequence) - window_size + 1),
                      seq_len(nchar(sequence) - window_size + 1) + window_size - 1)
  
  # Calculate properties for each window
  properties <- sapply(windows, function(x) {
    c(gc = calculate_gc_content(x),
      complexity = calculate_sequence_complexity(x),
      tm = calculate_melting_temp(x))
  })
  
  # Convert to long format for ggplot
  property_df <- as.data.frame(t(properties))
  property_df$position <- 1:nrow(property_df)
  property_df_long <- pivot_longer(property_df, 
                                 cols = c("gc", "complexity", "tm"),
                                 names_to = "property",
                                 values_to = "value")
  
  # Create heatmap using ggplot
  ggplot(property_df_long, aes(x = position, y = property, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = "DNA Properties",
         x = "Sequence Position",
         y = "Property",
         fill = "Value")
}

#' Create interactive 3D DNA structure visualization
#' @param sequence DNA sequence
#' @return Plotly object
create_3d_structure <- function(sequence) {
  # Calculate 3D coordinates
  n <- nchar(sequence)
  t <- seq(0, 10*pi, length.out = n)
  radius <- 2
  
  coords <- data.frame(
    x = radius * cos(t),
    y = radius * sin(t),
    z = seq(0, 10, length.out = n),
    nucleotide = strsplit(sequence, "")[[1]]
  )
  
  # Create backbone coordinates
  backbone1 <- data.frame(
    x = (radius - 0.2) * cos(t),
    y = (radius - 0.2) * sin(t),
    z = seq(0, 10, length.out = n)
  )
  
  backbone2 <- data.frame(
    x = (radius + 0.2) * cos(t),
    y = (radius + 0.2) * sin(t),
    z = seq(0, 10, length.out = n)
  )
  
  # Create 3D plot
  plot_ly() %>%
    add_trace(data = coords,
              x = ~x, y = ~y, z = ~z,
              color = ~nucleotide,
              colors = c("A" = "#FF9999", "T" = "#99FF99",
                        "G" = "#9999FF", "C" = "#FFFF99"),
              type = "scatter3d",
              mode = "markers",
              marker = list(size = 3)) %>%
    add_trace(data = backbone1,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d",
              mode = "lines",
              line = list(color = 'grey', width = 2),
              showlegend = FALSE) %>%
    add_trace(data = backbone2,
              x = ~x, y = ~y, z = ~z,
              type = "scatter3d",
              mode = "lines",
              line = list(color = 'grey', width = 2),
              showlegend = FALSE) %>%
    layout(title = "3D DNA Structure")
}

#' Create sequence logo using ggplot
#' @param sequences List of DNA sequences
#' @return ggplot object
create_sequence_logo <- function(sequences) {
  if(length(sequences) == 1) {
    # If single sequence, create frequency plot
    seq_chars <- strsplit(sequences, "")[[1]]
    freq_data <- data.frame(
      position = 1:length(seq_chars),
      nucleotide = seq_chars
    )
    
    ggplot(freq_data, aes(x = position, fill = nucleotide)) +
      geom_bar(stat = "count", position = "stack") +
      scale_fill_manual(values = c(
        "A" = "#FF9999", "T" = "#99FF99",
        "G" = "#9999FF", "C" = "#FFFF99"
      )) +
      theme_minimal() +
      labs(title = "Sequence Composition",
           x = "Position",
           y = "Count")
  } else {
    # Multiple sequences not implemented
    ggplot() +
      annotate("text", x = 0, y = 0, 
               label = "Multiple sequence logo requires\nalignment functionality") +
      theme_void()
  }
}

#' Create circular plot for DNA sequence
#' @param sequence DNA sequence
#' @return ggplot object
create_circular_plot <- function(sequence) {
  # Create data for circular layout
  bases <- strsplit(sequence, "")[[1]]
  n <- length(bases)
  angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
  radius <- 10
  
  plot_data <- data.frame(
    x = radius * cos(angles),
    y = radius * sin(angles),
    base = bases,
    position = 1:n
  )
  
  # Add connecting lines data
  lines_data <- data.frame(
    x = radius * c(cos(angles), cos(angles[-1]), NA),
    y = radius * c(sin(angles), sin(angles[-1]), NA),
    base = rep(bases, each = 2)
  )
  
  # Create plot
  ggplot() +
    geom_path(data = lines_data, aes(x = x, y = y, color = base), size = 0.5) +
    geom_point(data = plot_data, aes(x = x, y = y, color = base), size = 3) +
    scale_color_manual(values = c(
      "A" = "#FF9999",
      "T" = "#99FF99",
      "G" = "#9999FF",
      "C" = "#FFFF99"
    )) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "right") +
    labs(title = "Circular DNA Visualization", color = "Base")
}

#' Create secondary structure heatmap
#' @param pair_matrix Base pairing probability matrix
#' @return ggplot object
create_structure_heatmap <- function(pair_matrix) {
  # Convert matrix to long format
  n <- nrow(pair_matrix)
  df <- expand.grid(x = 1:n, y = 1:n)
  df$value <- as.vector(pair_matrix)
  
  # Create heatmap
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = "Base Pairing Probability",
         x = "Sequence Position",
         y = "Sequence Position")
}

#' Create codon usage plot
#' @param codon_data Codon analysis results
#' @return ggplot object
create_codon_plot <- function(codon_data) {
  df <- data.frame(
    codon = names(codon_data$frequencies),
    frequency = as.numeric(codon_data$frequencies)
  )
  
  ggplot(df, aes(x = reorder(codon, -frequency), y = frequency)) +
    geom_bar(stat = "identity", aes(fill = frequency)) +
    scale_fill_viridis_c() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Codon Usage Frequencies",
         x = "Codon",
         y = "Frequency")
}

#' Create flexibility plot
#' @param flex_data Flexibility analysis results
#' @return ggplot object
create_flexibility_plot <- function(flex_data) {
  ggplot(flex_data, aes(x = position, y = flexibility)) +
    geom_line() +
    geom_smooth(method = "loess", span = 0.1) +
    theme_minimal() +
    labs(title = "DNA Flexibility Profile",
         x = "Sequence Position",
         y = "Flexibility Score")
}

#' Create comprehensive property plot
#' @param sequence DNA sequence
#' @return plotly object
create_property_plot <- function(sequence) {
  # Calculate properties in sliding windows
  window_size <- 10
  n <- nchar(sequence) - window_size + 1
  
  properties <- data.frame(
    position = 1:n,
    gc = numeric(n),
    flex = numeric(n),
    complexity = numeric(n)
  )
  
  for(i in 1:n) {
    subseq <- substr(sequence, i, i + window_size - 1)
    properties$gc[i] <- calculate_gc_content(subseq)
    properties$complexity[i] <- calculate_sequence_complexity(subseq)
  }
  
  plot_ly() %>%
    add_trace(data = properties,
              x = ~position,
              y = ~gc,
              name = "GC Content",
              type = "scatter",
              mode = "lines") %>%
    add_trace(data = properties,
              x = ~position,
              y = ~complexity * 100,
              name = "Complexity",
              type = "scatter",
              mode = "lines") %>%
    layout(title = "DNA Properties",
           xaxis = list(title = "Position"),
           yaxis = list(title = "Value"))
}