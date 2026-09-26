# Calculate Accuracy Score

A simple Python implementation of the **accuracy score** for classification models using NumPy.

## Problem

Write a function that takes two 1D NumPy arrays:

- `y_true`: true labels
- `y_pred`: predicted labels

and returns the classification accuracy as a `float`.

## Formula

$$
\text{Accuracy} = \frac{\text{Number of Correct Predictions}}{\text{Total Number of Predictions}}
$$



## Explanation

The comparison

```python
y_true == y_pred
```

returns:

```text
[ True  True  False  True  True  True ]
```

There are **5 correct predictions** out of **6 total predictions**.

So:

$$
\text{Accuracy} = \frac{5}{6} = 0.8333
$$

or approximately:

$$
83.33\%
$$

## Requirements

- Python 3
- NumPy

Install NumPy with:

```bash
pip install numpy
```
