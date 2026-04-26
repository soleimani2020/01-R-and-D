# AdaBoost From Scratch in Python

## Overview

This repository contains a **from-scratch implementation of AdaBoost (Adaptive Boosting)** using **Python** and **NumPy**, built without external machine learning frameworks.

The project demonstrates how boosting combines multiple weak learners (decision stumps) into a strong classifier by iteratively adjusting sample weights and minimizing classification error.

---

## Key Features

- Pure Python + NumPy implementation
- Decision stump weak learners
- Feature, threshold, and polarity optimization
- Weighted error minimization
- Alpha weight calculation for each weak learner
- Iterative boosting rounds
- Educational and transparent codebase

---

## How AdaBoost Works

For each boosting iteration:

1. Initialize all sample weights equally
2. Search for the best decision stump:
   - Feature
   - Threshold
   - Polarity
3. Compute weighted classification error
4. Calculate learner weight:

```math
\alpha = \frac{1}{2} \ln\left(\frac{1 - error}{error}\right)
