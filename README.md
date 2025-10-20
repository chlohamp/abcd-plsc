# ABCD Sociocultural PLSC Analysis

This repository contains R scripts that used partial least squares correlation (PLSC) ([Krishnan et al., 2011](https://pubmed.ncbi.nlm.nih.gov/20656037/); [McIntosh and Lobaugh, 2004](https://pubmed.ncbi.nlm.nih.gov/15501095/_) to identify latent, multivariate patterns of associations between the adolescent sociocultural environment (at the individual, family, and community levels) and resting state functional connectivity (rsFC). rsFC was chosen as a neurobiological metric of interest given (i) prior work linking resting state fluctuations and cultural variables ([Constante et al., 2023](https://pubmed.ncbi.nlm.nih.gov/35429195/); [Luo et al., 2022](https://pubmed.ncbi.nlm.nih.gov/34160613/); [Meca et al., 2023](https://www.sciencedirect.com/science/article/pii/S2667174323000137)), (ii) the known functional relevance of rsFC to self-referential and affiliative processing ([Pintos Lobo et al., 2020](https://pubmed.ncbi.nlm.nih.gov/36436737/); [Raichle, 2015](https://pubmed.ncbi.nlm.nih.gov/25938726/_), and (iii) the flexibility of rsFC to examine the full range of brain systems in absence of external stimuli and/or task demands ([Smith et al., 2009](https://www.pnas.org/doi/10.1073/pnas.0905267106_). 

Additionally, for each significant latent dimension derived from PLSC, we conducted a regression analysis to determine how the identified sociocultural brain patterns may relate to behavioral risk/resilience factors, which will include measures of both physical health (e.g., physical inactivity, weight/obesity, sleep) and mental health (e.g., internalizing/externalizing symptoms).

## Python Environment Requirements

For the Python analysis components (including brain surface plotting and visualization), you will need:

- **Python 3.9.6**
- **Gradec**: Install first using:
  ```bash
  pip install git+https://github.com/JulioAPeraza/gradec.git
  ```

After installing Gradec, you can install the other required Python packages using the standard pip install commands for neuroimaging libraries (nilearn, surfplot, nibabel, etc.).

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

