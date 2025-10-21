# ABCD Sociocultural PLSC Analysis

This repository contains R scripts and Python notebooks that use partial least squares correlation (PLSC) ([Krishnan et al., 2011](https://pubmed.ncbi.nlm.nih.gov/20656037/); [McIntosh and Lobaugh, 2004](https://pubmed.ncbi.nlm.nih.gov/15501095/)) to identify latent, multivariate patterns of associations between the adolescent sociocultural environment (at the individual, family, and community levels) and resting state functional connectivity (rsFC). rsFC was chosen as a neurobiological metric of interest given (i) prior work linking resting state fluctuations and cultural variables ([Constante et al., 2023](https://pubmed.ncbi.nlm.nih.gov/35429195/); [Luo et al., 2022](https://pubmed.ncbi.nlm.nih.gov/34160613/); [Meca et al., 2023](https://www.sciencedirect.com/science/article/pii/S2667174323000137)), (ii) the known functional relevance of rsFC to self-referential and affiliative processing ([Pintos Lobo et al., 2020](https://pubmed.ncbi.nlm.nih.gov/36436737/); [Raichle, 2015](https://pubmed.ncbi.nlm.nih.gov/25938726/)), and (iii) the flexibility of rsFC to examine the full range of brain systems in absence of external stimuli and/or task demands ([Smith et al., 2009](https://www.pnas.org/doi/10.1073/pnas.0905267106)). 

Additionally, for each significant latent dimension derived from PLSC, we conducted a regression analysis to determine how the identified sociocultural brain patterns may relate to behavioral risk/resilience factors, which will include measures of both physical health (e.g., physical inactivity, weight/obesity, sleep) and mental health (e.g., internalizing/externalizing symptoms).

## Repository Structure

### Analysis Scripts

#### **1. clean-vars.ipynb**
**Purpose:** Data preprocessing and variable name standardization  
**Description:** This Jupyter notebook loads raw ABCD data files and standardizes variable names for clarity and brevity. It performs two main operations:
- **Sociocultural variables**: Renames measures to short codes (e.g., `meim_ss_exp` → `EIE-Y`, `via_ss_hc` → `HA-Y`)
- **RSFC variables**: Converts network connectivity column names from ABCD format to simplified network pairs (e.g., `rsfmri_c_ngd_cgc_ngd_dt` → `CON-DN`)

**Outputs:**
- `derivatives/none-reduced/clean-socio.csv` - Cleaned sociocultural data
- `derivatives/none-reduced/clean-rsfc.csv` - Cleaned RSFC data

#### **2. clean-plsc.Rmd**
**Purpose:** Main PLSC analysis pipeline  
**Description:** This R Markdown document performs the complete PLSC analysis workflow:
- **Residualization**: Removes effects of covariates (age, sex, scanner site, motion) from both RSFC and sociocultural data using multiple linear regression
- **PLSC Analysis**: Runs tepPLS to identify latent dimensions linking sociocultural measures to brain connectivity patterns
- **Permutation Testing**: Tests significance of dimensions using 10,000 permutations via `perm4PLSC`
- **Bootstrap Analysis**: Calculates bootstrap ratios to assess variable stability using `Boot4PLSC`
- **Variance Calculation**: Computes percent variance explained by each component
- **Visualization**: Creates scree plots, correlation plots, and bootstrap ratio bar plots

**Key Outputs:**
- Scree plot showing significant dimensions
- Correlation matrix heatmap
- Bootstrap ratio plots for each significant dimension
- CSV files with loadings and bootstrap ratios for both data blocks
- Latent variable scores (lx, ly)

**Required R Packages:**
- `TExPosition`, `PTCA4CATA`, `data4PCCAR` - PLSC functionality
- `tidyverse`, `dplyr` - Data manipulation
- `ggplot2`, `corrplot`, `cowplot` - Visualization

#### **3. plsc-figures.ipynb**
**Purpose:** Brain visualization and surface plotting  
**Description:** This Python notebook creates publication-ready brain visualizations of significant PLSC results:
- **Custom Colormaps**: Defines network-specific color schemes for each Gordon network
- **Individual Networks**: Creates surface plots for each significant network (Dimensions 1 & 3)
- **Combined Networks**: Generates overlapping network visualizations showing spatial relationships between network pairs
- **Composite Figures**: Assembles multi-panel figures for dimensions 1 and 3

**Key Features:**
- Uses `surfplot` and `neuromaps` for high-quality surface rendering on fsLR 32k surface
- Transforms Gordon parcellation from volumetric to surface space
- Creates both individual network plots and network pair combinations
- Generates publication-ready figures at 300 DPI

**Outputs:**
- Individual network surface plots (e.g., `CinguloOpercular_surface_plot.png`)
- Combined network plots (e.g., `CinguloOpercular_Default_combined_surface_plot.png`)
- Dimension-specific composite figures (`networks-dim1.png`, `networks-dim3.png`)
- All outputs saved to `derivatives/none-reduced/figures/`

#### **4. plsc-reg-corr.R**
**Purpose:** Regression analysis of PLSC network connectivity and health outcomes  
**Description:** This R script examines relationships between specific network connectivity pairs (identified from PLSC) and physical/mental health variables using linear mixed-effects models (LME). It tests associations while controlling for demographic covariates and accounting for nested data structure (participants within families within sites).

**Network Pairs Analyzed:**
- CON-DN (Cingulo-Opercular - Default Network)
- DN-DAN (Default - Dorsal Attention Network)
- DN-DN (Default - Default Network)
- DN-VN (Default - Visual Network)
- VN-VN (Visual - Visual Network)

**Health Variables Tested:**
- Sleep measures (sleep duration, sleep corrected)
- Respiratory health (wheeze, cough, diagnosis, bronchitis)
- Blood pressure (systolic, diastolic)
- Physical activity
- Behavioral problems (internalizing, externalizing symptoms)

**Statistical Model:**
```r
rsfc ~ health_variable + covariates + (1|site_id_l/rel_family_id)
```

**Key Features:**
- Controls for age, sex, parental demographics, income, scanner manufacturer, and head motion
- Accounts for nested structure with random effects for site and family
- Scales all continuous predictors for standardized coefficients
- Loops through all network pairs and health variables systematically

#### **5. plsc-reg-boot.R**
**Purpose:** Regression analysis of PLSC latent dimension scores and health outcomes  
**Description:** Similar to `plsc-reg-corr.R`, but analyzes the **latent variable scores** (composite brain-behavior patterns) from PLSC dimensions rather than individual network pairs. Tests how the overall brain-behavior patterns (Dimension 1 and Dimension 3) relate to health outcomes.

**Dimensions Analyzed:**
- **Dimension 1**: Primary latent pattern from PLSC
- **Dimension 3**: Secondary latent pattern from PLSC

**Statistical Model:**
```r
dimension_score ~ health_variable + covariates + (1|site_id_l/rel_family_id)
```

**Additional Features:**
- Creates visualizations for significant associations
- Uses `ggplot2` for publication-quality plots
- Focuses on significant variables (e.g., wheeze) for detailed plotting
- Provides holistic view of how brain-behavior patterns relate to health

**Comparison with plsc-reg-corr.R:**
- `plsc-reg-corr.R`: Tests specific **network connectivity pairs**
- `plsc-reg-boot.R`: Tests **overall latent dimension scores** (composite patterns)

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

1. **Data Preparation** (`clean-vars.ipynb`)
   - Load raw ABCD data
   - Standardize variable names
   - Save cleaned datasets

2. **PLSC Analysis** (`clean-plsc.Rmd`)
   - Residualize covariates
   - Run PLSC with permutation testing
   - Calculate bootstrap ratios
   - Generate statistical visualizations
   - Export latent variable scores and loadings

3. **Brain Visualization** (`plsc-figures.ipynb`)
   - Create custom network colormaps
   - Generate individual network surface plots
   - Create network pair combinations
   - Assemble multi-panel figures

4. **Health Outcomes Analysis** (R scripts)
   - **Network-level** (`plsc-reg-corr.R`): Test specific network connectivity pairs against health variables
   - **Pattern-level** (`plsc-reg-boot.R`): Test overall PLSC dimension scores against health variables
   - Control for demographic and technical covariates
   - Account for nested data structure (participants within families within sites)

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

