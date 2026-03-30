# ABCD Sociocultural PLSC Analysis

This repository contains R scripts and Python notebooks that use partial least squares correlation (PLSC) ([Krishnan et al., 2011](https://pubmed.ncbi.nlm.nih.gov/20656037/); [McIntosh and Lobaugh, 2004](https://pubmed.ncbi.nlm.nih.gov/15501095/)) to identify latent, multivariate patterns of associations between the adolescent sociocultural environment (at the individual, family, and community levels) and resting state functional connectivity (rsFC). rsFC was chosen as a neurobiological metric of interest given (i) prior work linking resting state fluctuations and cultural variables ([Constante et al., 2023](https://pubmed.ncbi.nlm.nih.gov/35429195/); [Luo et al., 2022](https://pubmed.ncbi.nlm.nih.gov/34160613/); [Meca et al., 2023](https://www.sciencedirect.com/science/article/pii/S2667174323000137)), (ii) the known functional relevance of rsFC to self-referential and affiliative processing ([Pintos Lobo et al., 2020](https://pubmed.ncbi.nlm.nih.gov/36436737/); [Raichle, 2015](https://pubmed.ncbi.nlm.nih.gov/25938726/)), and (iii) the flexibility of rsFC to examine the full range of brain systems in absence of external stimuli and/or task demands ([Smith et al., 2009](https://www.pnas.org/doi/10.1073/pnas.0905267106)). 

Additionally, for each significant latent dimension derived from PLSC, we conducted a regression analysis to determine how the identified sociocultural brain patterns may relate to behavioral risk/resilience factors, which will include measures of both physical health (e.g., physical inactivity, weight/obesity, sleep) and mental health (e.g., internalizing/externalizing symptoms).

## Repository Structure

### Analysis Scripts (in order of execution)

#### **0a-plsc-df.ipynb** - Extract ABCD Variables
**Purpose:** Extract raw ABCD Study data and organize into initial datasets  
**Description:** This Python notebook loads raw ABCD data files and extracts relevant variables for the analysis:
- **Sociocultural variables**: Ethnicity, identity, values, perceived discrimination, community factors from youth and parent-reported measures
- **Demographics**: Age, sex, parental education, income, ethnicity, family structure  
- **Linked external data (LED)**: Child Opportunity Index, neighborhood disadvantage, community statistics
- **Physical & mental health**: Sleep, respiratory function, blood pressure, physical activity, internalizing/externalizing symptoms
- **RSFC data**: Resting-state functional connectivity matrices from fMRI at 4-year follow-up
- **Imaging metadata**: Scanner manufacturer, MRI information, motion parameters

**Outputs:**
- `derivatives/none-reduced-no-motion/sociocult.csv` - Raw sociocultural variables
- `derivatives/none-reduced-no-motion/rsfc.csv` - Raw RSFC connectivity measures
- `derivatives/none-reduced-no-motion/covariate.csv` - Demographic covariates

#### **0c-vars-clean.ipynb** - Rename and Standardize Variables
**Purpose:** Clean variable names and create abbreviated naming convention  
**Description:** This Python notebook standardizes variable names for clarity and brevity:
- **Sociocultural variables**: Renames measures to short codes (e.g., `meim_ss_exp` → `EIE-Y` for ethnic identity exploration youth, `via_ss_hc` → `HA-Y` for heritage affirmation youth)
- **RSFC variables**: Converts network connectivity column names from ABCD format to simplified network pairs (e.g., `rsfmri_c_ngd_dt_ngd_dt` → `DN-DN` for Default Network-Default Network)
- **Creates cleaned dataframes** with consistent naming for downstream PLSC analysis

**Outputs:**
- `derivatives/none-reduced-no-motion/clean-socio.csv` - Cleaned sociocultural data with abbreviated names
- `derivatives/none-reduced-no-motion/clean-rsfc.csv` - Cleaned RSFC data with network pair names

#### **0d-demographics.ipynb** - Generate Sample Demographics
**Purpose:** Summarize demographic characteristics of the analyzed sample  
**Description:** Creates demographic summary tables and descriptive statistics:
- Age distribution (interview age of participants)
- Sex distribution
- Racial/ethnic composition  
- Parental education and income levels
- Site representation
- Family structure information

**Outputs:**
- Summary statistics tables
- Demographic breakdowns by key variables

#### **1a-plsc.Rmd** - (Deprecated - Do Not Use)
**Purpose:** Earlier version of PLSC pipeline (superseded by 1b)  
**Note:** This script is no longer used in the current analysis workflow. Use **1b-plsc-clean.Rmd** instead.

#### **1b-plsc-clean.Rmd** - Main PLSC Analysis Pipeline
**Purpose:** Complete Partial Least Squares Correlation (PLSC) analysis  
**Description:** This R Markdown document performs the complete PLSC analysis workflow to identify latent multivariate patterns linking sociocultural environment to brain connectivity:
- **Residualization**: Removes effects of covariates (age, sex, scanner site, motion, parental factors, income) from both RSFC and sociocultural data using multiple linear regression models
- **PLSC Analysis**: Runs `tepPLS()` to identify latent dimensions linking sociocultural measures to brain connectivity patterns via Singular Value Decomposition of the correlation matrix
- **Permutation Testing**: Tests significance of dimensions using 10,000 permutations via `perm4PLSC()` to determine which dimensions are significant at p < 0.01
- **Bootstrap Analysis**: Calculates bootstrap ratios to assess stability and reliability of variables using `Boot4PLSC()` with 10,000 iterations
- **Variance Calculation**: Computes percent variance explained by each significant component
- **Visualization**: Creates scree plots, correlation heatmaps, and bootstrap ratio bar plots

**Key Statistical Outputs:**
- Eigenvalues and p-values for each dimension
- Percent variance explained by each component
- Loadings (projections) for both data blocks (RSFC and sociocultural)
- Bootstrap ratios indicating variable stability
- Latent variable scores for significant dimensions

**Outputs:**
- Scree plot showing significant vs. non-significant dimensions
- Correlation matrix heatmap between data blocks
- Bootstrap ratio bar plots for each significant dimension (Dim 1, 2, 3, 4)
- CSV files with loadings: `rni-brain_P-components.csv`, `rni-sociocult_Q-components.csv`
- Latent scores: `rniXrsfc_lx-base.csv`, `rniXrsfc_ly-base.csv`
- Percent variance: `percent_variance_explained.csv`

**Required R Packages:**
- `TExPosition`, `PTCA4CATA`, `data4PCCAR` - PLSC functionality
- `tidyverse`, `dplyr` - Data manipulation
- `ggplot2`, `corrplot`, `cowplot` - Visualization
- `kableExtra` - Table formatting

#### **1c-plsc-figures.ipynb** - Brain Surface Visualizations
**Purpose:** Create publication-ready brain surface visualizations  
**Description:** This Python notebook generates high-quality brain surface plots showing which networks and brain regions form the significant PLSC dimensions:
- **Custom Network Colormaps**: Defines network-specific color schemes for 13 Gordon networks based on typical color assignments
- **Individual Network Plots**: Creates surface plots for each significant network pair (e.g., Default Network, Visual Network, etc.)
- **Combined Network Visualizations**: Generates overlapping surface plots showing spatial relationships of network pairs identified in PLSC
- **Composite Multi-panel Figures**: Assembles final publication-ready figures for Dimensions 1 and 3
- **High-resolution Output**: Renders figures at 300 DPI for publication quality

**Key Features:**
- Uses `surfplot` and `neuromaps` for high-quality surface rendering on fsLR 32k surface template
- Transforms Gordon parcellation from volumetric to surface space
- Creates consistent color schemes across all figures
- Generates both individual network plots and network pair combinations

**Outputs:**
- Individual network surface plots (e.g., `Default_surface_plot.png`, `Visual_surface_plot.png`)
- Combined network pair plots (e.g., `Default_Visual_combined_surface_plot.png`)
- Dimension-specific composite figures (`networks-dim1.png`, `networks-dim3.png`)
- All outputs saved to `derivatives/none-reduced-no-motion/figures/`

#### **2a-regression-df.ipynb** - Prepare Regression Analysis Data
**Purpose:** Create analysis-ready dataframes for regression studies  
**Description:** This Python notebook prepares data for testing whether PLSC-identified brain patterns predict physical health outcomes:
- **Filters subjects**: Retains only participants present in physical health dataset
- **Creates regression datasets**: For each significant network pair and latent dimension, creates a CSV with:
  - The specific RSFC measure (for network pairs) or latent score (for dimensions)
  - All demographic covariates (age, sex, site, family structure, etc.)
  - All physical health outcome variables (sleep, respiratory, blood pressure, activity, behavioral measures)
- **Creates separate dim1 and dim3 folders**: Organizes data by dimension for parallel analysis
- **Saves subject-filtered versions** of all datasets with consistent subject IDs across files

**Network Pairs (Dimension 1) analyzed:**
- DN-DN (Default Network - Default Network)
- DN-VN (Default Network - Visual Network)
- VN-VN (Visual Network - Visual Network)

**Network Pairs (Dimension 3) analyzed:**
- DN-SMN (Default Network - Sensorimotor Mouth Network)

**Outputs:**
- `derivatives/none-reduced-no-motion/regression/dim1/phyhealth_*.csv` - Data for regression models
- `derivatives/none-reduced-no-motion/regression/dim3/phyhealth_*.csv` - Data for regression models
- Filtered versions of: `rsfc-reg.csv`, `covariate-reg.csv`, `latent-reg.csv`, `phyhealth-reg-mod.csv`

#### **2b-plsc-reg-corr.R** - Regression: Network Connectivity × Health Outcomes
**Purpose:** Test if specific RSFC network pairs predict physical health outcomes  
**Description:** This R script examines relationships between specific network connectivity pairs (identified from PLSC as significant) and physical/mental health variables using linear mixed-effects models (LME):
- **Model specification**: Tests association between each RSFC network pair and health variable while controlling for demographic covariates and accounting for nested data structure
- **Random effects**: Includes random intercepts for site and family (nested within site) to account for data hierarchy
- **Fixed effects**: Adjusts for age, sex, parental demographics, income, scanner manufacturer, and MRI technical factors
- **Categorical variables**: Properly handles categorical predictors (sex, origin, manufacturer) using sum-to-zero contrasts for Type III ANOVA
- **Loops through both dimensions**: Automatically processes all dim1 networks and dim3 networks specified in scripts

**Statistical Models:**
```r
RSFC_measure ~ Health_Variable + Age + Sex + Parental_Age + Parental_Gender + 
               Parental_Education + Partner_Education + Income + Origin + 
               Site + Scanner + (1|Site/Family)
```

**Summary Outputs:**
- `plsc-reg-corr-results-dim1.csv` - Results for dimension 1 networks
- `plsc-reg-corr-results-dim3.csv` - Results for dimension 3 networks  
- `plsc-reg-corr-results-all.csv` - Combined results
- Per-model coefficient tables saved for significant associations

**Key Outputs:**
- p-values for each network-health variable association
- Significance flags (p < 0.05, p < 0.01)
- Model coefficients and standard errors
- Diagnostic plots for model assumptions

#### **2b-plsc-reg-boot.R** - Regression: Latent Dimensions × Health Outcomes
**Purpose:** Test if PLSC latent dimension scores predict physical health outcomes  
**Description:** Similar to `2b-plsc-reg-corr.R`, but analyzes the **overall latent variable scores** from PLSC dimensions rather than individual network pairs. Tests how the composite brain-behavior patterns relate to health outcomes:
- **Analyzes latent dimension scores**: Uses V1 (Dimension 1) and V3 (Dimension 3) latent scores as predictors
- **Same LME structure**: Controls for demographics and accounts for nested data (site/family)
- **Composite patterns**: Provides holistic view comparing overall brain-behavior associations to health
- **Parallel structure**: Organized to mirror `2b-plsc-reg-corr.R` for comparison

**Statistical Models:**
```r
Latent_Score ~ Health_Variable + Age + Sex + Parental_Age + Parental_Gender + 
               Parental_Education + Partner_Education + Income + Origin + 
               Site + Scanner + (1|Site/Family)
```

**Outputs:**
- `plsc-reg-boot-results-dim1.csv` - Results for dimension 1 latent scores
- `plsc-reg-boot-results-dim3.csv` - Results for dimension 3 latent scores
- `plsc-reg-boot-results-all.csv` - Combined results
- Diagnostic visualizations for significant effects

**Comparison of 2b scripts:**
| Aspect | plsc-reg-corr.R | plsc-reg-boot.R |
|--------|-----------------|-----------------|
| Predictor | Individual network pairs | Overall latent dimension scores |
| Level of analysis | Network-specific | Pattern-level (composite) |
| Use case | Which specific networks matter? | Does the overall brain pattern matter? |
| Output granularity | Detailed by network pair | Holistic by dimension |

#### **2c-regression-figures.ipynb** - Visualize Health Outcome Associations
**Purpose:** Create publication-ready plots of significant regression findings  
**Description:** This Python notebook visualizes significant associations between PLSC-identified brain patterns and physical health outcomes:
- **Scatter plots with regression lines**: Shows bivariate relationships between RSFC/latent scores and health outcomes
- **Conditional plots by group**: Displays associations separated by categorical covariates (e.g., by sex, site)
- **Effect size visualization**: Highlightssignificance and magnitude of associations
- **Multi-panel figures**: Assembles figures comparing dimensions or networks
- **Publication-ready formatting**: High resolution, consistent styling, manuscript-ready annotations

**Outputs:**
- Individual association plots for significant network-health pairs
- Dimension comparison plots
- Multi-panel composite figures
- All figures saved to `derivatives/none-reduced-no-motion/regression/figures/`

## Environment Setup

### Python Environment

For the Python analysis components (including brain surface plotting and visualization), you will need:

- **Python 3.9.6**
- **Gradec**: Install first using:
  ```bash
  pip install git+https://github.com/JulioAPeraza/gradec.git
  ```

After installing Gradec, install the required neuroimaging packages:
```bash
pip install nilearn surfplot neuromaps nibabel matplotlib seaborn pillow
```

### R Environment

For the PLSC analysis in R:

- **R version 4.5.1** (or 4.3.x/4.4.x)
- **Required packages** (install from CRAN):
  ```r
  install.packages(c("tidyverse", "dplyr", "ggplot2", "corrplot", 
                     "cowplot", "readr", "gridExtra", "psych", 
                     "car", "kableExtra", "ExPosition"))
  ```

- **GitHub packages** (install using devtools):
  ```r
  library(devtools)
  install_github("HerveAbdi/TExPosition")
  install_github("HerveAbdi/PTCA4CATA")
  install_github("HerveAbdi/data4PCCAR")
  ```

## Workflow

The analysis follows a sequential, hierarchical pipeline from raw data extraction through statistical validation:

### Stage 1: Data Preparation (0x scripts)
1. **0a-plsc-df.ipynb** - Extract raw ABCD variables (sociocultural, demographics, health, imaging)
2. **0c-vars-clean.ipynb** - Rename and standardize variable names for analysis
3. **0d-demographics.ipynb** - Generate demographic summary tables

*Inputs:* Raw ABCD data files from `dset/` directory  
*Outputs:* Clean, standardized datasets ready for statistical analysis

### Stage 2: PLSC Analysis (1x scripts)  
1. **1a-plsc.Rmd** - (Deprecated - skip this script)
2. **1b-plsc-clean.Rmd** - Perform complete PLSC analysis pipeline:
   - Residualize covariates from both data blocks
   - Identify latent dimensions linking sociocultural environment to brain connectivity
   - Determine statistical significance via permutation testing (10,000 iterations)
   - Assess variable stability via bootstrap analysis (10,000 iterations)
   - Generate statistical visualizations (scree plots, heatmaps, bootstrap ratio plots)
3. **1c-plsc-figures.ipynb** - Create publication-ready brain surface visualizations

*Inputs:* Clean datasets from Stage 1  
*Outputs:* Significant dimensions, loadings, latent scores, brain surface plots, percent variance explained

### Stage 3: Regression Analysis (2x scripts)
1. **2a-regression-df.ipynb** - Prepare analysis dataframes:
   - Filter subjects to those with physical health data
   - Create separate datasets for each significant network pair and latent dimension
   - Merge with all covariates and health outcome variables
   - Organize into dim1 and dim3 subdirectories

2. **2b-plsc-reg-corr.R** - Test regression relationships between network pairs and health outcomes:
   - Linear mixed effects models accounting for nested data structure (site/family)
   - Controls for age, sex, demographics, scanner manufacturer, imaging parameters
   - Generates p-values, effect sizes, and model diagnostics
   - Outputs separate results for dim1 and dim3

3. **2b-plsc-reg-boot.R** - Test regression relationships between latent dimension scores and health outcomes:
   - Uses same LME framework but with overall dimension scores as predictors
   - Provides holistic view of brain-behavior-health associations
   - Parallel structure to 2b-plsc-reg-corr.R for direct comparison

4. **2c-regression-figures.ipynb** - Visualize significant health outcome associations:
   - Create scatter plots with regression lines
   - Show conditional associations by demographic groups
   - Assemble multi-panel publication-ready figures

*Inputs:* PLSC results from Stage 2, health data from Stage 1  
*Outputs:* Statistical test results, significance tables, effect size visualizations

## Variables

Variables that were included in the PLSC analysis:

| Sociocultural Measures          | Gordon Networks                  |
|----------------------------------|----------------------------------|
| MEIM-R youth                      | Auditory Network                 |
| VIA youth                          | Cingulo-opercular Network        |
| MACV youth                          | Cingulo-parietal Network         |
| MEIM-R caregiver                    | Default Network                  |
| VIA caregiver                        | Dorsal Attention Network         |
| MACV caregiver                        | Fronto-parietal Network          |
| Perceived discrimination              | Retrosplenial Temporal Network   |
| Community cohesion                     | Sensorimotor Hand Network        |
| LED: Child opportunity index           | Sensorimotor Mouth Network       |
| LED: NaNDA disadvantage index          | Salience Network                 |
| LED: Getis-Ord Gi* Statistic           | Ventral Attention Network        |
|                                          | Visual Network                   |

