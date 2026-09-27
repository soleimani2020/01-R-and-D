# Linear Regression Using the Normal Equation

## Overview

Linear regression finds coefficients that best describe the relationship between input features `X` and a target variable `y`.

Using the **Normal Equation**, the coefficients can be calculated directly without gradient descent:

\[
\theta = (X^T X)^{-1} X^T y
\]

where:

- `X` = feature matrix
- `y` = target vector
- `X.T` = transpose of `X`
- `θ` = regression coefficients

---
