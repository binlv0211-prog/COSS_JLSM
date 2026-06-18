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
*   **Real Data 3 (`real_data3.R`, `real_data3_functions.R`):** Analysis of the **Teenage Friends and Lifestyle Study**.

---

## Environment and Usage

*   **Language:** All scripts are written in **R**.
*   **Requirements:** To execute the code, ensure you have a functional R environment and install the necessary dependencies (e.g., `MASS`, `pgdraw`, or other packages invoked in the scripts).
