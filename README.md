# GenomeScope - DNA Sequence Analyzer

A comprehensive web-based tool for DNA sequence analysis built with R Shiny. GenomeScope provides researchers and students with an intuitive interface for analyzing nucleotide sequences, calculating molecular properties, and visualizing genomic data.

## Features

**Core Analysis**
- Nucleotide composition and frequency analysis
- GC content calculation with visual representation
- Melting temperature prediction
- Sequence validation and error checking
- Complementary sequence generation

**Advanced Analysis**
- GC skew analysis across sequence regions
- Sequence complexity calculations
- Pattern recognition and repeat identification
- Codon usage statistics
- Secondary structure predictions

**Visualizations**
- Interactive nucleotide frequency charts
- 3D DNA structure modeling
- Circular genome representations
- Property distribution heatmaps
- Statistical summary tables

## Screenshots

### Main Dashboard
Main dashboard interface featuring intuitive navigation and modern design
![Dashboard Overview](screenshots/Dashboard%20About%20Page.png)

### Sequence Input Interface
User-friendly sequence input with validation and example loading
![Sequence Input](screenshots/Sequence%20Input%20&%20Validation.png)

### Basic Analysis Results
Comprehensive nucleotide composition analysis with interactive visualizations
![Nucleotide Analysis](screenshots/Nucleotide%20Composition%20Analysis.png)

### Detailed Composition Analysis
![Detailed Analysis](screenshots/Nucleotide%20Composition%20Analysis%202.png)

### Advanced Analysis Tools
Advanced sequence analysis including GC skew, complexity patterns, and distribution plots
![Advanced Analysis](screenshots/Advanced%20Analysis%20.png)

### 3D Structure Visualization
Interactive 3D DNA double helix visualization with customizable viewing angles
![3D Structure](screenshots/3D%20Structure%20Analysis.png)

### Statistical Reports
Detailed sequence statistics, pattern recognition, and molecular properties
![Statistics](screenshots/Detailed%20Statistics%20&%20Reports.png)

## Installation

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)

### Required Packages
```r
# Install CRAN packages
install.packages(c(
    "shiny", "shinydashboard", "stringr", "ggplot2", 
    "plotly", "DT", "tidyr", "dplyr", "viridis"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "GenomicRanges", "seqinr"))
```

### Running the Application
```bash
# Clone the repository
git clone https://github.com/[YOUR_USERNAME]/dna-analyser-r-shiny.git
cd dna-analyser-r-shiny

# Open R/RStudio and run
shiny::runApp()
```

## Usage

1. **Input**: Enter or paste DNA sequences in the input field
2. **Validate**: The system automatically validates sequence format
3. **Analyze**: Click the analyze button to process the sequence
4. **Explore**: Navigate through different analysis tabs to view results
5. **Export**: Download charts and statistical reports

## Technical Details

**Built With**
- R Shiny for the web interface
- Biostrings for sequence manipulation
- ggplot2 and Plotly for visualizations
- DT for interactive data tables

**Analysis Methods**
- Sliding window analysis for regional properties
- Nearest-neighbor calculations for melting temperature
- Pattern matching algorithms for repeat detection
- Statistical methods for sequence complexity

## Project Structure
```
dna-analyser-r-shiny/
├── app.R                           # Main application file
├── R/
│   ├── dna_functions.R            # Core analysis functions
│   └── visualization_helpers.R    # Plotting utilities
├── www/
│   └── about.md                   # Documentation
├── tests/
│   └── testthat/                  # Unit tests
└── screenshots/                   # Application screenshots
```

## Contributing

Contributions are welcome. Please feel free to submit issues or pull requests to improve the functionality or fix bugs.


## Contact

**Developer**: Winnie Cook
**Email**: winniecook19@gmail.com 
**GitHub**: [@winniecook](https://github.com/winniecook)
