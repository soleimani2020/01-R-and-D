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

Calculate R-squared for Regression Analysis

Formula

$$
R^2 = 1 - \frac{SS_{\text{fit}}}{SS_{\text{mean}}}
$$

where:

$$
SS_{\text{fit}} = \sum_{i=1}^{n}(y_i-\hat{y}_i)^2
$$

and:

$$
SS_{\text{mean}} = \sum_{i=1}^{n}(y_i-\bar{y})^2
$$

So:

SS(fit) measures the squared error of the regression model.

SS(mean) measures the squared error we would get if we predicted the mean of y for every observation.

Therefore:

$$
R^2 = 1 - \frac{\text{Error of fitted model}}{\text{Error of mean model}}
$$
