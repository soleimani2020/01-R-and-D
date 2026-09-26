# Calculate Root Mean Square Error (RMSE)

A simple Python implementation of the **Root Mean Square Error (RMSE)** metric using NumPy.

## Problem

Write a function:

```python
rmse(y_true, y_pred)
```

that calculates the RMSE between:

- `y_true`: actual values
- `y_pred`: predicted values

The function should return the RMSE rounded to **three decimal places**.

It should also handle:

- Mismatched array shapes
- Empty arrays
- Invalid input types

## Formula

$$
RMSE = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2}
$$

where:

- $y_i$ = actual value
- $\hat{y}_i$ = predicted value
- $n$ = number of observations

RMSE first squares the prediction errors, calculates their mean, and then takes the square root.


