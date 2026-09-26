# Implement F-Score Calculation for Binary Classification

A simple Python implementation of the **F-Score** for binary classification using NumPy.

## Problem

Write a function:

```python
f_score(y_true, y_pred, beta)
```

where:

- `y_true`: true binary labels
- `y_pred`: predicted binary labels
- `beta`: controls the importance of Precision versus Recall

The function should return the F-Score rounded to **three decimal places**.

## Formula

$$
F_\beta =
(1+\beta^2)
\frac{\text{Precision}\times\text{Recall}}
{\beta^2\times\text{Precision}+\text{Recall}}
$$

When `beta = 1`, this becomes the **F1-Score**:

$$
F_1 =
2
\frac{\text{Precision}\times\text{Recall}}
{\text{Precision}+\text{Recall}}
$$

where:

$$
\text{Precision} = \frac{TP}{TP+FP}
$$

and:

$$
\text{Recall} = \frac{TP}{TP+FN}
$$

## Interpretation of Beta

- `beta = 1`: Precision and Recall have equal importance
- `beta > 1`: Recall is given more importance
- `beta < 1`: Precision is given more importance

