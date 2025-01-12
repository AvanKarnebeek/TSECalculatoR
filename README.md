# TSECalculatoR

## Overview

For the **Bachelor End Thesis Evaluating the T Cell-to-Stroma Enrichment (TSE) Score as a Transcriptomic Predictive Biomarker for Immunotherapy Response**, the R package `TSECalculatoR` was created. This package allows users to extract the TSE score, TSE category, and additional biomarkers from a count matrix, making it a valuable tool for transcriptomic data analysis.

## Features

The `TSECalculatoR` package consists of three core functions, all requiring the same input: a **count matrix** with genes as rows and patient IDs as columns.

### 1. `calculate_tsescore(countmatrix)`
This function uses the Gene Set Variation Analysis (GSVA) method to:
- Calculate the **TSE score**.
- Extract the **global T cell signature**.
- Extract the **global stromal cell signature**.

### 2. `calculate_tseprofexh(countmatrix)`
This function also uses GSVA to:
- Calculate the **TSE score**.
- Extract the **global T cell and stromal cell signatures**.
- Calculate a **proliferation score**.
- Calculate a **T cell exhaustion score**.

### 3. `TSE_classify(countmatrix)`
This function:
- Classifies patients into one of three TSE categories (**Negative**, **Neutral**, or **Positive**) based on centroids.
- Returns the **correlation values** to each centroid.
- Provides the **p-value** of the correlation.
- This function is adapted from the code developed by **Rijnders et al.**, available at [https://github.com/ANakauma/TSEscore_ICIs](https://github.com/ANakauma/TSEscore_ICIs).

## Installation

You can install the `TSECalculatoR` package directly from GitHub using the following command:

```R
# Install devtools if not already installed
install.packages("devtools")

# Install TSECalculatoR from GitHub
devtools::install_github("AvanKarnebeek/TSECalculatoR")
```

Replace `AvanKarnebeek` with your GitHub username.

## Usage

Here’s an example of how to use the functions in `TSECalculatoR`:

```R
# Load the package
library(TSECalculatoR)

# Example: Load your count matrix
data <- read.csv("path_to_your_count_matrix.csv", row.names = 1)

# Calculate TSE score and signatures
results_tse <- calculate_tsescore(data)

# Calculate TSE score, proliferation, and exhaustion scores
results_prof_exh <- calculate_tseprofexh(data)

# Classify patients into TSE categories
classification_results <- TSE_classify(data)

# View the results
head(results_tse)
head(results_prof_exh)
head(classification_results)
```

## Input Requirements
- A **count matrix** with:
  - Genes as rows.
  - Patient IDs as columns.

## Citation
If you use this package, please cite the code developed by **Rijnders et al.** available at [https://github.com/ANakauma/TSEscore_ICIs](https://github.com/ANakauma/TSEscore_ICIs).

## License
This package is released under the [MIT License](LICENSE).

---

For any questions or issues, please contact Alexine van Karnebeek at a.vankarnebeek@erasmusmc.nl . 
