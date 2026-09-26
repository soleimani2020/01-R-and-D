# Implement Recall Metric in Binary Classification

A simple Python implementation of the **Recall** metric for binary classification.

## Problem

Write a function:

```python
recall(y_true, y_pred)
```

that calculates recall from:

- `y_true`: true binary labels
- `y_pred`: predicted binary labels

The function should return a `float`.

If `TP + FN = 0`, return `0.0` to avoid division by zero.

## Formula

$$
\text{Recall} = \frac{TP}{TP + FN}
$$

where:

- **TP (True Positives):** predicted `1` and actually `1`
- **FN (False Negatives):** predicted `0` but actually `1`

