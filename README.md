# Description

This repository contains the R code required to reproduce the numerical findings presented in the following paper:

> **Lv, B., Tang, Y., & Zhang, S. (2025).** *A Cumulative Ordered Spike-and-Slab Prior for Adaptive Dimension Selection in Joint Latent Space Models.*
---

## Core Functionality 

The repository provides a comprehensive framework for joint latent space modeling, including data generation, model estimation, and performance evaluation. All core functions used in the implementation are organized in the `functions/` folder.

### 1. Data Simulation and Evaluation
*   `data_process.R`: A utility script for generating synthetic datasets and calculating evaluation metrics (e.g., recovery accuracy, predictive performance).

### 2. Gibbs Sampling Algorithms
A suite of scripts implementing Gibbs samplers for various data configurations:
*   **Marginal Models:**
    *   `network_only.R`: Network structure only.
    *   `Y_only_normal.R` / `binary.R` / `ordinal.R` / `mix.R`: Node attributes only (Continuous, Binary, Ordinal, or Mixed-type).
*   **Joint Models (Network + Node Attributes):**
    *   `network_normal.R`: Network data with continuous attributes.
    *   `network_binary.R`: Network data with binary attributes.
    *   `network_ordinal.R`: Network data with ordinal attributes.
    *   `network_mix.R`: Network data with mixed-type attributes (continuous, binary, and ordinal).

### 3. Cross-Validation and Model Selection
Scripts for dimension selection and prior sensitivity analysis using normal priors and Kullback-Leibler (KL) divergence:
*   `network_normal_kl.R`: Cross-validation for continuous attributes in networks.
*   `network_binary_kl.R`: Cross-validation for binary attributes in networks.

### 4. Missing Data Imputation
Specialized functions for handling missingness in binary node attributes:
*   **Latent Space Approaches:** `network_miss.R` and `network_binary_miss.R` leverage the latent space structure for imputation.
*   **Baseline Comparison:** `mice.R` implements the Multiple Imputation by Chained Equations (MICE) algorithm.

---

## Replication Scripts

High-level scripts are provided to reproduce the core experiments and case studies:

### Simulation Studies
*   **Sample Size (`Simulation1.R`):** Investigates model performance across varying sample sizes.
*   **Network Density (`Simulation2.R`):** Explores the impact of network density within a six-node simulation setup.
*   **Hyperparameter Sensitivity (`Simulation_sensitive.R`):** Evaluates the influence of hyperparameter choices ($a, a_\theta, b_\theta$).
*   **Attribute Types (`Simulation_ordinal.R`, `Simulation_mix.R`):** Validates the model's performance on ordinal and mixed-type node attributes.

### Empirical Applications
*   **Real Data 1 (`real_data1.R`):** Analysis of the **French Financial Elites** dataset.
*   **Real Data 2 (`real_data2.R`):** Application to the **Facebook Social Network** dataset.
*   **Real Data 3 (`real_data3/`):** Analysis of the **Teenage Friends and Lifestyle Study**. The analysis pipeline includes: `01_functions.R` (core functions), `02_fit_data.R` (model fitting), `03_dim_select.R` (dimension selection), `04_metric_cal.R` (metric calculation), and `05_plot.ipynb` (visualization in Jupyter Notebook).

---

## Data Access and Preprocessing

The empirical applications in this repository are based on three publicly available benchmark datasets. Data sources and preprocessing steps are summarized below.

### 1. French Financial Elites

The French Financial Elites dataset was originally collected by Kadushin (1995) to study social connections among members of the French financial elite during the final years of France's Socialist government.

For this study, we use the processed version distributed with the **`jlsm`** R package (Wang, 2021), which has also been analyzed in Wang et al. (2022) and Wang et al. (2023).

**Data access:**

The processed dataset can be obtained from the archived source of the **`jlsm`** package:

* CRAN Archive: https://cran.r-project.org/src/contrib/Archive/jlsm/

The package contains the dataset `french`, which includes the friendship network and associated node attributes used in our analysis.

**Preprocessing:**

No additional preprocessing was performed beyond the version distributed with the `jlsm` package.


---

### 2. Facebook Ego Networks

The Facebook Social Network dataset consists of ego-networks collected from Facebook users and is publicly available through the Stanford Network Analysis Project (SNAP).

**Data access:**

https://snap.stanford.edu/data/

**Preprocessing performed in this repository:**

* We analyze 8 of the 10 ego-networks.
* Two networks are excluded due to extremely small size or highly sparse attribute vectors.
* Binary attributes with extremely low or high prevalence are removed within each network, following Zhang et al. (2022).
* The remaining anonymized binary attributes are used for missing-data imputation experiments.

---

### 3. Teenage Friends and Lifestyle Study (Glasgow Study)

The Teenage Friends and Lifestyle Study is a longitudinal study of adolescent friendship networks and health-related behaviors conducted in Glasgow between 1995 and 1997.

The dataset is distributed through the SIENA project.

**Data access:**

https://www.stats.ox.ac.uk/~snijders/siena/Glasgow_data.htm

**Preprocessing performed in this repository:**

* We use the older cohort at the first measurement wave.
* Only the 129 students observed at all three waves are retained.
* Friendship nominations are converted to an undirected binary network:

  * `A[i,j] = 1` if either student nominated the other.
* The original distinction between “best friends” and “just friends” is collapsed into a single friendship indicator.
* Node attributes include ordinal and mixed-type variables describing substance use, leisure activities, and musical preferences.

---


## Environment and Usage

*   **Language:** All scripts are written in **R**, with the exception of visualization components (e.g., `05_plot.ipynb` in `real_data3/`) which are implemented in **Python**.
*   **Requirements:** To execute the code, ensure you have functional **R** and **Python** environments. Install the necessary R dependencies (e.g., `MASS`, `pgdraw`, or other packages invoked in the scripts) and Python dependencies (e.g., `matplotlib`, `pandas`, `numpy`) for running the Jupyter notebooks.
