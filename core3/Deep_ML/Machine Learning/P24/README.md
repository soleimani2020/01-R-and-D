# Calculate R-squared for Regression Analysis

A simple Python implementation of the **R-squared (R²)** metric for regression using NumPy.

## Problem

Write a function:

```python
r_squared(y_true, y_pred)
```

that calculates the **R-squared** value from:

- `y_true`: actual target values
- `y_pred`: predicted target values

The result should be rounded to **three decimal places**.

## Formula

$$
R^2 = 1 - \frac{SS_{res}}{SS_{tot}}
$$

where:

$$
SS_{res} = \sum_{i=1}^{n}(y_i - \hat{y}_i)^2
$$

and:

$$
SS_{tot} = \sum_{i=1}^{n}(y_i - \bar{y})^2
$$

- $SS_{res}$ = residual sum of squares
- $SS_{tot}$ = total sum of squares
- $y_i$ = actual value
- $\hat{y}_i$ = predicted value
- $\bar{y}$ = mean of the actual values



```bash
pip install numpy
```
