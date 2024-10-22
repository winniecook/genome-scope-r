# DNA Sequence Analyzer - R Shiny App

## Overview
An advanced DNA sequence analysis tool built with R Shiny, providing interactive visualization and comprehensive analysis of DNA sequences. This application offers both basic and advanced DNA sequence analysis features through an intuitive web interface.
RShiny version of my previous dna sequence analyser app built in python/javascript

## Features

### Basic Analysis
- Nucleotide frequency analysis with interactive visualization
- GC content calculation and pie chart visualization
- Codon usage statistics and distribution
- Basic sequence statistics

### Advanced Analysis
- Nucleotide distribution across sequence
- GC skew analysis
- Sequence complexity analysis
- Pattern recognition and repeat identification

### Structural Analysis
- Interactive 3D DNA structure visualization
- Circular DNA visualization
- DNA property analysis including:
  - Melting temperature
  - Flexibility
  - GC content distribution
  - Molecular properties

## Installation

### Prerequisites
- R (>= 4.0.0)
- RStudio (recommended for development)

### Required R Packages
```R
# Install required CRAN packages
install.packages(c(
    "shiny",
    "shinydashboard",
    "stringr",
    "ggplot2",
    "plotly",
    "DT",
    "tidyr",
    "dplyr",
    "viridis",
    "markdown"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges"))
```

### Setup
1. Clone the repository:
```bash
git clone https://github.com/yourusername/dna-analyzer-r-shiny.git
cd dna-analyzer-r-shiny
```

2. Open RStudio and set the working directory to the project folder

3. Run the application:
```R
shiny::runApp()
```

## Usage

1. Input Tab:
   - Enter or paste your DNA sequence
   - Use the "Load Example" button for a sample sequence
   - Click "Analyze" to process the sequence

2. Basic Analysis Tab:
   - View nucleotide frequencies
   - Examine GC content distribution
   - Analyze codon usage patterns

3. Advanced Analysis Tab:
   - Explore nucleotide distribution patterns
   - Analyze GC skew
   - Investigate sequence complexity

4. Structure Tab:
   - Interact with 3D DNA structure
   - View circular DNA representation
   - Examine DNA physical properties

## Project Structure
```
dna-analyzer-r-shiny/
├── R/
│   ├── dna_functions.R      # Core analysis functions
│   └── visualization_helpers.R  # Visualization functions
├── www/
│   └── about.md            # About page content
├── app.R                   # Main application file
└── README.md              # Project documentation
```

## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
- Built with R Shiny framework
- Uses Biostrings for DNA sequence analysis
- Visualization powered by ggplot2 and plotly
