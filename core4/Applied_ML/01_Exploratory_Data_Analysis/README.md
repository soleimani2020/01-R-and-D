# Project 01: Kernel Ridge Regression (KRR) for SANS Inverse Problems

**Category:** Scientific Machine Learning  
**Methods:** Kernel Ridge Regression, inverse modeling, regression, model evaluation  
**Application:** Small-Angle Neutron Scattering (SANS)

---

## Overview

This project investigates the application of Kernel Ridge Regression (KRR) to inverse problems in Small-Angle Neutron Scattering (SANS).

The objective is to develop a machine learning workflow that maps scattering intensity profiles to underlying physical model parameters. Starting with simulated scattering data, the project explores model training, hyperparameter optimization, prediction accuracy, and robustness to noise.

## Objectives

- Generate or prepare synthetic SANS scattering datasets with known physical parameters.
- Develop a KRR regression model using scikit-learn.
- Optimize kernel and regularization hyperparameters.
- Evaluate parameter recovery using appropriate regression metrics.
- Investigate the effects of noise and measurement uncertainty.
- Examine parameter identifiability and the limitations of data-driven inverse modeling.

## Methodology

1. **Data Generation:** Generate synthetic SANS scattering curves using a well-defined physical model.
2. **Data Preprocessing:** Prepare scattering intensity profiles and corresponding physical parameters.
3. **Baseline Model:** Implement Kernel Ridge Regression using scikit-learn.
4. **Hyperparameter Optimization:** Use cross-validation to optimize the kernel and regularization parameters.
5. **Model Evaluation:** Compare predicted parameters with their known values using MAE, RMSE, and other appropriate metrics.
6. **Robustness Analysis:** Investigate the effects of noise and instrumental resolution.
7. **Scientific Interpretation:** Examine parameter identifiability and assess the limitations of the inverse model.

## Physical Model

The initial implementation will focus on spherical particles.

Potential parameters to investigate include:

- Particle radius
- Polydispersity
- Scattering intensity scale
- Background intensity

The initial model will focus on a small subset of these parameters to establish a reliable baseline.

## Technologies

- Python
- NumPy
- pandas
- Matplotlib
- scikit-learn
- Jupyter Notebook

## Project Structure

```text
01_KRR_SANS_Inverse_Problems/
├── README.md
├── main.py
└── project.ipynb
