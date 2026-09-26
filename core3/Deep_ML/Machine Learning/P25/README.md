# Calculate Mean Absolute Error (MAE)

A simple Python implementation of the **Mean Absolute Error (MAE)** metric using NumPy.

## Problem

Write a function:

```python
mean_absolute_error(y_true, y_pred)
```

that calculates the average absolute difference between the actual and predicted values.

- `y_true`: actual values
- `y_pred`: predicted values

The function should return the MAE as a `float`.

## Formula

$$
MAE = \frac{1}{n}\sum_{i=1}^{n}|y_i - \hat{y}_i|
$$

where:

- $y_i$ = actual value
- $\hat{y}_i$ = predicted value
- $n$ = total number of observations

