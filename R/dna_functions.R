#' DNA Analysis Functions
#' This file contains core functions for DNA sequence analysis

library(Biostrings)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)
library(seqinr)
library(GenomicRanges)

#' Validate DNA sequence
#' @param sequence Character string containing DNA sequence
#' @return Boolean indicating if sequence is valid DNA
is_valid_dna <- function(sequence) {
  # Convert to uppercase and remove whitespace
  sequence <- toupper(str_trim(sequence))
  # Check if sequence contains only valid DNA characters
  return(all(str_split(sequence, "")[[1]] %in% c("A", "T", "C", "G")))
}

#' Calculate nucleotide frequencies
#' @param sequence Character string containing DNA sequence
#' @return Named vector of nucleotide counts
count_nucleotides <- function(sequence) {
  sequence <- toupper(sequence)
  counts <- c(
    A = str_count(sequence, "A"),
    T = str_count(sequence, "T"),
    G = str_count(sequence, "G"),
    C = str_count(sequence, "C")
  )
  return(counts)
}

#' Calculate GC content
#' @param sequence Character string containing DNA sequence
#' @return Numeric value representing GC percentage
calculate_gc_content <- function(sequence) {
  counts <- count_nucleotides(sequence)
  gc_content <- (counts["G"] + counts["C"]) / sum(counts) * 100
  return(round(gc_content, 2))
}

#' Generate complementary sequence
#' @param sequence Character string containing DNA sequence
#' @return Character string of complementary sequence
get_complementary_sequence <- function(sequence) {
  # Convert to DNAString object
  dna <- DNAString(sequence)
  # Get complement
  comp <- as.character(complement(dna))
  return(comp)
}

#' Plot nucleotide frequencies
#' @param counts Named vector of nucleotide counts
#' @return ggplot object
plot_nucleotide_frequencies <- function(counts) {
  # Create plot data
  plot_data <- data.frame(
    Nucleotide = names(counts),
    Count = as.numeric(counts)
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = Nucleotide, y = Count, fill = Nucleotide)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    scale_fill_manual(values = c(
      "A" = "#FF9999",
      "T" = "#99FF99",
      "G" = "#9999FF",
      "C" = "#FFFF99"
    )) +
    labs(title = "Nucleotide Frequencies")
  
  return(p)
}

#' Plot GC content
#' @param gc_content Numeric GC percentage
#' @return ggplot object
plot_gc_content <- function(gc_content) {
  # Create data
  plot_data <- data.frame(
    Content = c("GC", "AT"),
    Percentage = c(gc_content, 100 - gc_content)
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = "", y = Percentage, fill = Content)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = c(
      "GC" = "#9999FF",
      "AT" = "#FF9999"
    )) +
    labs(title = "GC Content")
  
  return(p)
}

#' Calculate DNA sequence properties
#' @param sequence Character string containing DNA sequence
#' @return List of physical properties
calculate_dna_properties <- function(sequence) {
  # Calculate various properties
  properties <- list(
    length = nchar(sequence),
    molecular_weight = sum(sapply(strsplit(sequence, "")[[1]], function(x) {
      switch(x,
             "A" = 331.2,
             "T" = 322.2,
             "G" = 347.2,
             "C" = 307.2)
    })),
    melting_temp = calculate_melting_temp(sequence),
    complexity = calculate_sequence_complexity(sequence)
  )
  
  return(properties)
}

#' Calculate melting temperature
#' @param sequence DNA sequence
#' @return Numeric melting temperature
calculate_melting_temp <- function(sequence) {
  # Basic nearest-neighbor method
  counts <- count_nucleotides(sequence)
  gc_content <- calculate_gc_content(sequence)
  
  # Simplified formula for sequences < 14 bp
  if(nchar(sequence) < 14) {
    tm <- (counts['A'] + counts['T']) * 2 + (counts['G'] + counts['C']) * 4
  } else {
    # Wallace formula for longer sequences
    tm <- 64.9 + 41 * (gc_content / 100 - 0.41)
  }
  
  return(round(tm, 2))
}

#' Calculate sequence complexity
#' @param sequence DNA sequence
#' @return Numeric complexity score
calculate_sequence_complexity <- function(sequence) {
  # Calculate linguistic complexity
  k_max <- min(nchar(sequence), 10)
  observed_kmers <- sapply(1:k_max, function(k) {
    length(unique(substring(sequence,
                          seq_len(nchar(sequence) - k + 1),
                          seq_len(nchar(sequence) - k + 1) + k - 1)))
  })
  
  possible_kmers <- sapply(1:k_max, function(k) {
    min(4^k, nchar(sequence) - k + 1)
  })
  
  complexity <- sum(observed_kmers) / sum(possible_kmers)
  return(round(complexity, 4))
}

#' Find repeating patterns
#' @param sequence DNA sequence
#' @param min_length Minimum pattern length
#' @return Data frame of repeating patterns
find_repeating_patterns <- function(sequence, min_length = 3) {
  patterns <- data.frame(
    pattern = character(),
    start = integer(),
    length = integer(),
    count = integer(),
    stringsAsFactors = FALSE
  )
  
  for(i in min_length:min(20, nchar(sequence))) {
    # Get all possible patterns of length i
    all_patterns <- substring(sequence,
                            seq_len(nchar(sequence) - i + 1),
                            seq_len(nchar(sequence) - i + 1) + i - 1)
    
    # Find repeating patterns
    repeat_patterns <- table(all_patterns)
    repeat_patterns <- repeat_patterns[repeat_patterns > 1]
    
    if(length(repeat_patterns) > 0) {
      for(pattern in names(repeat_patterns)) {
        positions <- as.integer(str_locate_all(sequence, pattern)[[1]][,1])
        patterns <- rbind(patterns, data.frame(
          pattern = pattern,
          start = I(list(positions)),
          length = nchar(pattern),
          count = as.integer(repeat_patterns[pattern]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(patterns)
}

#' Predict DNA structural properties
#' @param sequence DNA sequence
#' @return List of structural properties
predict_dna_structure <- function(sequence) {
  # Calculate base stacking energies
  dinucleotides <- substring(sequence,
                           seq_len(nchar(sequence) - 1),
                           seq_len(nchar(sequence) - 1) + 1)
  
  # Stacking energy values (kcal/mol)
  stacking_energies <- list(
    "AA" = -5.37, "AT" = -6.57, "AG" = -6.78, "AC" = -6.57,
    "TA" = -3.82, "TT" = -5.37, "TG" = -6.57, "TC" = -5.37,
    "GA" = -6.78, "GT" = -6.57, "GG" = -8.26, "GC" = -9.61,
    "CA" = -6.57, "CT" = -6.78, "CG" = -9.61, "CC" = -8.26
  )
  
  total_energy <- sum(sapply(dinucleotides, function(x) {
    stacking_energies[[x]] %||% 0
  }))
  
  return(list(
    stacking_energy = total_energy,
    bendability = 0, # Placeholder
    stability = -total_energy/nchar(sequence)
  ))
}

#' Analyze DNA sequence (main function)
#' @param sequence Character string containing DNA sequence
#' @return List of analysis results
analyze_dna <- function(sequence) {
  # Basic analyses
  counts <- count_nucleotides(sequence)
  gc <- calculate_gc_content(sequence)
  comp <- get_complementary_sequence(sequence)
  
  return(list(
    nucleotide_counts = counts,
    gc_content = gc,
    complementary_sequence = comp
  ))
}

#' Create advanced plots for DNA sequence
#' @param sequence DNA sequence
#' @return List of plot objects
create_advanced_plots <- function(sequence) {
  # Nucleotide distribution along sequence
  sliding_window <- function(seq, window = 10) {
    windows <- substring(seq,
                        seq_len(nchar(seq) - window + 1),
                        seq_len(nchar(seq) - window + 1) + window - 1)
    
    counts <- t(sapply(windows, function(x) {
      c(A = str_count(x, "A"),
        T = str_count(x, "T"),
        G = str_count(x, "G"),
        C = str_count(x, "C"))
    }))
    
    return(as.data.frame(counts))
  }
  
  window_data <- sliding_window(sequence)
  window_data$position <- seq_len(nrow(window_data))
  
  # Create plots
  plots <- list(
    distribution = ggplot(pivot_longer(window_data, -position),
                         aes(x = position, y = value, color = name)) +
      geom_line() +
      theme_minimal() +
      labs(title = "Nucleotide Distribution",
           x = "Position",
           y = "Count",
           color = "Nucleotide"),
    
    gc_skew = ggplot(window_data,
                     aes(x = position,
                         y = (G - C)/(G + C))) +
      geom_line() +
      theme_minimal() +
      labs(title = "GC Skew",
           x = "Position",
           y = "GC Skew")
  )
  
  return(plots)
}



#' Predict DNA secondary structure
#' @param sequence DNA sequence
#' @return List of structural predictions
predict_secondary_structure <- function(sequence) {
  # Calculate base pairing probability
  n <- nchar(sequence)
  bases <- strsplit(sequence, "")[[1]]
  
  # Simple base pairing rules
  pair_scores <- matrix(0, n, n)
  for(i in 1:(n-4)) {
    for(j in (i+4):n) {
      # Check complementary bases
      if((bases[i] == "A" && bases[j] == "T") ||
         (bases[i] == "T" && bases[j] == "A") ||
         (bases[i] == "G" && bases[j] == "C") ||
         (bases[i] == "C" && bases[j] == "G")) {
        pair_scores[i,j] <- 1
      }
    }
  }
  
  return(list(
    pair_matrix = pair_scores,
    stability_score = sum(pair_scores)/n,
    potential_hairpins = find_hairpins(sequence)
  ))
}

#' Find potential hairpin structures
#' @param sequence DNA sequence
#' @return Data frame of potential hairpins
find_hairpins <- function(sequence) {
  hairpins <- data.frame(
    start = integer(),
    end = integer(),
    loop_size = integer(),
    stem_length = integer(),
    sequence = character(),
    stringsAsFactors = FALSE
  )
  
  n <- nchar(sequence)
  bases <- strsplit(sequence, "")[[1]]
  
  for(i in 1:(n-10)) {
    for(loop_size in 3:7) {
      for(stem_length in 3:6) {
        if(i + 2*stem_length + loop_size > n) next
        
        stem1 <- bases[i:(i+stem_length-1)]
        stem2 <- rev(bases[(i+stem_length+loop_size):(i+2*stem_length+loop_size-1)])
        
        # Check complementarity
        if(all(stem1 == complement(stem2))) {
          hairpins <- rbind(hairpins, data.frame(
            start = i,
            end = i + 2*stem_length + loop_size - 1,
            loop_size = loop_size,
            stem_length = stem_length,
            sequence = substr(sequence, i, i + 2*stem_length + loop_size - 1),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(hairpins)
}

#' Analyze codon usage
#' @param sequence DNA sequence
#' @return List of codon statistics
analyze_codon_usage <- function(sequence) {
  # Ensure sequence length is divisible by 3
  n <- nchar(sequence)
  n <- floor(n/3) * 3
  sequence <- substr(sequence, 1, n)
  
  # Split into codons
  codons <- substring(sequence,
                     seq(1, n-2, by=3),
                     seq(3, n, by=3))
  
  # Count codons
  codon_counts <- table(codons)
  
  # Calculate frequencies
  frequencies <- prop.table(codon_counts)
  
  # Identify start and stop codons
  start_positions <- which(codons == "ATG")
  stop_positions <- which(codons %in% c("TAA", "TAG", "TGA"))
  
  return(list(
    codon_counts = codon_counts,
    frequencies = frequencies,
    start_codons = start_positions,
    stop_codons = stop_positions,
    gc_content_third = calculate_gc_third_position(codons)
  ))
}

#' Calculate GC content at third codon position
#' @param codons Vector of codons
#' @return Numeric GC percentage
calculate_gc_third_position <- function(codons) {
  third_positions <- substr(codons, 3, 3)
  gc_count <- sum(third_positions %in% c("G", "C"))
  return((gc_count/length(codons)) * 100)
}

#' Predict DNA flexibility
#' @param sequence DNA sequence
#' @return Data frame of flexibility scores
predict_dna_flexibility <- function(sequence) {
  # Dinucleotide flexibility parameters (arbitrary scale)
  flex_params <- list(
    "AA" = 0.6, "AT" = 0.9, "AG" = 0.7, "AC" = 0.6,
    "TA" = 0.9, "TT" = 0.6, "TG" = 0.8, "TC" = 0.7,
    "GA" = 0.7, "GT" = 0.8, "GG" = 0.5, "GC" = 0.6,
    "CA" = 0.6, "CT" = 0.7, "CG" = 0.8, "CC" = 0.5
  )
  
  # Calculate flexibility for each position
  dinucs <- substring(sequence,
                     seq_len(nchar(sequence)-1),
                     seq_len(nchar(sequence)-1) + 1)
  
  flex_scores <- sapply(dinucs, function(x) flex_params[[x]])
  
  return(data.frame(
    position = 1:(length(flex_scores)),
    flexibility = flex_scores
  ))
}