# Close_relationship_and_role_of_fractal_dimension_in-unperturbed-tumor-growth
Computational framework for analyzing unperturbed tumor growth through mKJMA, MBBG, and Richards models, focusing on fractal dimension, heterogeneous nucleation, dynamic transformations, cooperativity, and phase transitions

## Overview

This repository contains the computational framework, numerical simulations, and data-analysis tools developed to investigate the relationship among **fractal dimension, heterogeneous nucleation, dynamic transformations, cooperativity mechanisms, and unperturbed tumor growth**.

The methodology provides an integrated mathematical description of tumor growth kinetics based on three growth models:

* Modified Kolmogorov–Johnson–Mehl–Avrami (**mKJMA**) model
* Montijano–Bergues–Bory–Gompertz (**MBBG**) model
* Richards model

The computational framework establishes mathematical relationships among these models and their characteristic parameters. By combining normalized growth kinetics, local dynamical indicators, susceptibility analysis, and fractal geometry, the project investigates tumor growth as a **stepwise, multiscale, and intrinsically fractal process involving multiple cooperative phase transitions**.

---

## Scientific Motivation

Tumor growth is commonly described using phenomenological growth laws. However, individual growth models often emphasize different aspects of the same biological process, making it difficult to establish a unified interpretation of tumor growth kinetics.

This project investigates whether different mathematical growth models can be connected through their parameters and whether these relationships reveal common underlying mechanisms.

The central hypothesis is that unperturbed tumor growth is governed by the interplay of:

* fractal geometry;
* heterogeneous nucleation;
* cooperative phenomena;
* competitive dynamics;
* dynamic transformations;
* phase-transition-like behavior;
* soft and hard impingement mechanisms.

The framework combines global growth models with local dynamical indicators to identify changes that may remain hidden in conventional tumor growth curves.

---

## Mathematical Models

### Modified Kolmogorov–Johnson–Mehl–Avrami Model

The modified KJMA framework is used to describe tumor growth in terms of transformation kinetics, heterogeneous nucleation, growth mechanisms, and impingement effects.

The model provides a connection between tumor growth kinetics and concepts originally developed in the theory of nucleation and phase transformations.

The analysis includes the estimation of the **local Avrami exponent**, which acts as a dynamical indicator of changes in the underlying growth mechanism.

---

### Montijano–Bergues–Bory–Gompertz Model

The MBBG model provides a generalized description of sigmoidal tumor growth dynamics.

In this repository, the model is normalized to facilitate direct comparison with the modified KJMA and Richards models. The normalization procedure enables the analysis of relationships among characteristic parameters and growth regimes.

---

### Richards Model

The Richards model provides a flexible generalized growth law capable of representing different asymmetric sigmoidal growth trajectories.

Its normalized formulation is used to investigate the correspondence between phenomenological tumor growth models and transformation-based kinetic descriptions.

---

## Model Correspondence

One of the main objectives of this project is to investigate the mathematical correspondence among the **mKJMA, MBBG, and Richards models**.

The computational analysis focuses on:

* comparison of normalized tumor growth curves;
* relationships among characteristic model parameters;
* correspondence between growth regimes;
* local changes in growth kinetics;
* fractal properties of tumor morphology;
* cooperative and competitive dynamics;
* transition points between dynamical regimes.

The results support a unified interpretation in which apparently different growth models describe related aspects of the same underlying multiscale tumor dynamics.

---

## Local Dynamical Indicators

### Local Avrami Exponent

The local Avrami exponent is used to detect temporal changes in tumor growth mechanisms.

Instead of assuming a single constant exponent throughout the complete growth process, the local formulation enables the identification of changes associated with:

* heterogeneous nucleation;
* growth dimensionality;
* cooperative transformations;
* competition among growing regions;
* impingement effects;
* transitions between avascular and vascular growth regimes.

---

### Susceptibility

A susceptibility-like order parameter is used to identify changes in the collective tumor growth dynamics.

The susceptibility analysis reveals different dynamical regimes throughout the growth trajectory and provides evidence of:

* cooperative heterogeneous nucleation during avascular growth;
* critical dynamical transformations;
* competitive behavior during vascular growth;
* transitions between cooperative and competitive regimes.

Together, the local Avrami exponent and susceptibility provide complementary tools for detecting regime changes during the complete tumor growth kinetics.

---

## Fractal Analysis

Tumor morphology is analyzed through the **contour fractal dimension**.

The computational tools included in this repository allow the analysis of changes in geometrical complexity throughout tumor growth. Discontinuities and variations in the fractal dimension are interpreted together with the local dynamical indicators.

The analysis suggests that tumor growth is intrinsically fractal and that morphological changes are associated with dynamical transitions in the underlying growth mechanisms.

The framework therefore connects:

* tumor morphology;
* growth kinetics;
* fractal geometry;
* phase dynamics;
* cooperative transformations.

---

## Main Results

The numerical simulations and data analysis support the following conclusions:

* The mKJMA, MBBG, and Richards models exhibit mathematical correspondence after appropriate normalization.
* Their characteristic parameters are interconnected and display relationships associated with the fractal nature of tumor growth.
* Complete unperturbed tumor growth kinetics is characterized by multiple dynamical regimes.
* Changes in contour fractal dimension are associated with transitions in growth dynamics.
* The local Avrami exponent detects temporal changes in transformation mechanisms.
* Susceptibility analysis reveals cooperative heterogeneous nucleation during avascular growth.
* Competitive dynamical regimes emerge during vascular growth.
* Tumor growth exhibits multiple cooperative phase transitions.
* Impingement mechanisms play an important role in the evolution of tumor growth kinetics.
* The complete growth process can be interpreted as stepwise, multiscale, and intrinsically fractal.

---

## Repository Structure

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
│
├── data/
│   ├── raw/
│   └── processed/
│
├── src/
│   ├── mKJMA/
│   ├── MBBG/
│   ├── Richards/
│   ├── local_avrami/
│   ├── susceptibility/
│   └── fractal_analysis/
│
├── simulations/
│   ├── growth_curves/
│   ├── model_comparison/
│   └── parameter_analysis/
│
├── figures/
│   ├── growth_kinetics/
│   ├── local_avrami/
│   ├── susceptibility/
│   └── fractal_dimension/
│
└── scripts/
    ├── reproduce_figures/
    └── statistical_analysis/
```

The exact directory structure may be adapted according to the programming language and organization of the simulation codes.

---

## Main Features

* Numerical implementation of the modified KJMA model
* Numerical implementation of the MBBG model
* Numerical implementation of the Richards model
* Normalization of tumor growth curves
* Parameter estimation and comparison
* Local Avrami exponent calculation
* Susceptibility analysis
* Detection of dynamical regime changes
* Contour fractal dimension analysis
* Identification of cooperative and competitive regimes
* Visualization of multiscale tumor growth kinetics
* Reproduction of manuscript figures

---

## Requirements

The numerical codes can be implemented and executed using the software environment associated with the repository.

Example requirements:

```text
MATLAB R2022a or later
Statistics and Machine Learning Toolbox
Optimization Toolbox
Image Processing Toolbox
```

If Python scripts are included, a typical environment may require:

```text
Python >= 3.10
NumPy
SciPy
Pandas
Matplotlib
scikit-image
```

Please adapt this section according to the final software dependencies of the repository.

---

## Usage

Clone the repository:

```bash
git clone <repository-url>
cd <repository-name>
```

For MATLAB implementations, add the source directories to the MATLAB path:

```matlab
addpath(genpath('src'));
```

Run the main analysis script:

```matlab
run('scripts/main_analysis.m');
```

To reproduce the figures associated with the manuscript:

```matlab
run('scripts/reproduce_figures/reproduce_all_figures.m');
```

The exact script names should be updated according to the final repository structure.

---

## Reproducing the Analysis

A typical computational workflow consists of the following steps:

1. Load experimental or simulated tumor growth data.
2. Normalize the tumor growth curves.
3. Fit the mKJMA, MBBG, and Richards models.
4. Estimate the corresponding model parameters.
5. Analyze relationships among parameters.
6. Calculate the local Avrami exponent.
7. Estimate the susceptibility order parameter.
8. Calculate the contour fractal dimension.
9. Identify transitions between dynamical regimes.
10. Compare avascular and vascular tumor growth phases.
11. Generate the figures and statistical analyses.

---

## Applications

The computational framework may be useful in the following research areas:

* Mathematical oncology
* Tumor growth modeling
* Cancer biophysics
* Nonlinear dynamical systems
* Fractal geometry
* Heterogeneous nucleation
* Phase transitions
* Cooperative phenomena
* Growth kinetics
* Cancer biomechanics
* Multiscale biological modeling
* Complex systems

---

## Data

The `data/` directory contains the datasets required to reproduce the analyses presented in the associated study.

When applicable, data should be organized into:

```text
data/
├── raw/          # Original datasets
└── processed/    # Normalized or transformed datasets
```

The original data should remain unchanged in the `raw/` directory. All transformations and preprocessing procedures should be performed through reproducible scripts.

---

## Citation

If you use the code, methodology, numerical procedures, or data-analysis tools contained in this repository, please cite the associated publication:

> **Close Relationship and Role of Fractal Dimension, Heterogeneous Nucleation, Dynamic Transformations, and Cooperativity Mechanisms in Unperturbed Tumor Growth**

BibTeX template:

```bibtex
@article{AuthorYearTumorGrowth,
  title   = {Close Relationship and Role of Fractal Dimension,
             Heterogeneous Nucleation, Dynamic Transformations
             and Cooperativity Mechanisms in Unperturbed Tumor Growth},
  author  = {Author(s)},
  journal = {Journal Name},
  year    = {Year},
  volume  = {Volume},
  pages   = {Pages},
  doi     = {DOI}
}
```

Please replace the bibliographic fields with the final publication information.

---

## License

This project is distributed under the license specified in the `LICENSE` file.

Please consult the license terms before reusing, modifying, or redistributing the source code.

---

## Contact

For questions about the mathematical models, numerical simulations, or scientific methodology, please contact the corresponding author of the associated publication.

---

## Keywords

`Tumor Growth` · `Mathematical Oncology` · `Fractal Dimension` · `mKJMA` · `MBBG` · `Richards Model` · `Heterogeneous Nucleation` · `Cooperativity` · `Phase Transitions` · `Avrami Exponent` · `Susceptibility` · `Cancer Biophysics` · `Multiscale Dynamics`
