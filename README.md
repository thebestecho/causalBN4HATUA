# Modelling Interventions to Combat Antibacterial Resistance in East Africa Using Causal Bayesian Networks

[![DOI](https://zenodo.org/badge/956904441.svg)](https://doi.org/10.5281/zenodo.21595803)

This repository provides the official source code, high-performance computing (HPC) deployment scripts, and analytical workflows for the manuscript:
> **Modelling interventions to combat antibacterial resistance in East Africa using causal Bayesian networks**  
> Published in *Communications Medicine*.

---

## 📌 Repository Architecture & Directory Structure

The repository is structured into three main sequential modules covering hyperparameter optimization, network structure learning, and causal intervention analysis:

```text
causalBN4HATUA/
├── simulation/
│   └── look_for_best_iss.R           # Hyperparameter optimization (ISS) via synthetic DAG simulation
├── BN_structure_learning/
│   ├── overall.R                     # Bootstrap Structural EM for BN structure learning & data imputation
│   └── run_structure_array.sh        # SLURM HPC Cluster array job script for structure learning
├── causal_BN_analysis/
│   ├── process_network_array.py      # CausalNex-based Bayesian Network fitting & Do-calculus intervention
│   └── run_causal_array.sh           # SLURM HPC Cluster array job script for counterfactual simulations
└── README.md                         # Project documentation
```

---

## 🛠️ Environment & Dependencies

### R Requirements (Modules 1 & 2)
Required R packages (**R >= 4.0.0**):
* [`bnlearn`](https://www.bnlearn.com/) - Structure learning with Structural EM (`structural.em`) and score optimization (`score`)
* [`boot`](https://cran.r-project.org/web/packages/boot/index.html) - Non-parametric bootstrap resampling
* [`MCMCpack`](https://cran.r-project.org/web/packages/MCMCpack/) - Dirichlet distributions for CPT sampling
* [`readr`](https://readr.tidyverse.org/), [`dplyr`](https://dplyr.tidyverse.org/), [`tidyr`](https://tidyr.tidyverse.org/), [`purrr`](https://purrr.tidyverse.org/), [`magrittr`](https://magrittr.tidyverse.org/) - Data handling and piping
* [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/html/00Index.html), [`doSNOW`](https://cran.r-project.org/web/packages/doSNOW/) - Multiprocessing support

Install R dependencies:
```R
install.packages(c("bnlearn", "boot", "MCMCpack", "readr", "dplyr", "tidyr", "purrr", "magrittr", "parallel", "doSNOW"))
```

### Python Requirements (Module 3)
Required Python packages (**Python >= 3.8**):
* [`causalnex`](https://github.com/quantumblacklabs/causalnex) - Structure models, Bayesian Networks, and Pearl's Do-calculus inference engine
* [`pandas`](https://pandas.pydata.org/) - Data manipulation and CSV exports

Install Python dependencies:
```bash
pip install causalnex pandas
```

---

## 🔬 Analytical Pipeline & Detailed Functions

### 1. `simulation/` — Hyperparameter Optimization
* **`look_for_best_iss.R`**: 
  * Generates ground-truth random Directed Acyclic Graphs (DAGs) using Ide & Cozman's `ic-dag` algorithm.
  * Parameterizes networks via Dirichlet distributions (`rdirichlet`) and samples synthetic datasets.
  * Evaluates BDe scores across an Imaginary Sample Size (`iss`) grid (1–100) to locate the optimal hyperparameter setting (determined as `iss = 18`).

---

### 2. `BN_structure_learning/` — Bootstrap Structure Learning with Structural EM
Handles missing data imputation and network structure learning across 1,000 bootstrap iterations on HPC clusters.

* **`overall.R`**:
  * **Bootstrap Resampling**: Takes array task IDs to randomly sample observations with replacement (`bootstrap_data`).
  * **Structural EM**: Executes `structural.em` using Hill-Climbing (`hc`) and BDe scoring (`iss = 18`, `maxp = 3`) constrained by a domain prior blacklist (`bl`).
  * **Outputs**: Saves learned DAG structures (`.RData`) to `./overall_structures/` and EM-imputed datasets (`.csv`) to `./overall_imputed/`.
* **`run_structure_array.sh`**:
  * SLURM array job script executing parallel instances of `overall.R`.

---

### 3. `causal_BN_analysis/` — Causal Inference & Counterfactual Interventions
Fits conditional probability distributions (CPDs) and simulates hypothetical exogenous public health interventions using Pearl's do-calculus.

* **`process_network_array.py`**:
  * **Network Reconstruction**: Reconstructs the DAG structure using `causalnex.structure.StructureModel` from arc list CSVs.
  * **CPD Parameter Fitting**: Fits node states and estimates CPDs using `BayesianEstimator` with K2 priors (`bayes_prior="K2"`).
  * **Stratified Baseline Rates**: Computes observational cross-tabulation for the target intervention variable vs. Multidrug Resistance (`MDR`).
  * **Do-Calculus Interventions (`InferenceEngine`)**: Applies `ie.do_intervention()` to simulate counterfactual scenarios where 100% of the population assumes a specific category state.
  * **Outputs**: Calculates query probabilities for `MDR` under baseline and counterfactual scenarios, saving results to `./overall_analysis_results/<variable>/results_<file_index>.csv`.

* **`run_causal_array.sh`**:
  * SLURM cluster array job script submitting 1,000 parallel evaluation jobs across 8 target intervention variables (`age`, `education`, `workingstatus`, `deprived_assets`, `overcrowded`, `toilet`, `protect_drinkingwater`, `aware_AMR`).

---

## 🚀 Reproduction Steps

### Step 1: Optimize Hyperparameters
```bash
cd simulation/
Rscript look_for_best_iss.R
```

### Step 2: Learn Network Structures (SLURM HPC)
```bash
cd BN_structure_learning/
sbatch run_structure_array.sh
```

### Step 3: Run Counterfactual Interventions (SLURM HPC)
```bash
cd causal_BN_analysis/
sbatch run_causal_array.sh
```

---

## 📄 License & Citation

This codebase is made available under the [MIT License](LICENSE).

If you use this repository in your research, please cite:
```bibtex
@article{Xuejia2026CausalBN,
  title={Modelling interventions to combat antibacterial resistance in East Africa using causal Bayesian networks},
  author={Xuejia et al.},
  journal={Communications Medicine},
  year={2026}
}
```
