# Implement Precision Metric

A simple Python implementation of the **Precision** metric for binary classification using NumPy.

## Problem

Write a Python function `precision` that calculates the precision metric given two NumPy arrays:

- `y_true`: true binary labels
- `y_pred`: predicted binary labels

Precision measures how many of the samples predicted as **positive** are actually positive.

## Formula

$$
\text{Precision} = \frac{TP}{TP + FP}
$$

where:

- **TP (True Positives):** predicted `1` and actually `1`
- **FP (False Positives):** predicted `1` but actually `0`


## Explanation

For the example:

```text
y_true = [1, 0, 1, 1, 0, 1]
y_pred = [1, 1, 0, 1, 0, 1]
```

The model predicted `1` four times.

Among these predictions:

```text
True Positives  = 3
False Positives = 1
```

Therefore:

$$
\text{Precision} = \frac{3}{3+1} = 0.75
$$

So the precision is:

```text
75%
```

## How the NumPy Conditions Work

True positives are found using:

```python
(y_true == 1) & (y_pred == 1)
```

False positives are found using:

```python
(y_true == 0) & (y_pred == 1)
```

`np.sum()` then counts the number of `True` values.


```bash
pip install numpy
```
