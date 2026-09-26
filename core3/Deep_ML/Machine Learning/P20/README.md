# Calculate Accuracy Score

A simple Python implementation of the **accuracy score** for classification models using NumPy.

## Problem

Write a function that takes two 1D NumPy arrays:

- `y_true`: true labels
- `y_pred`: predicted labels

and returns the classification accuracy as a `float`.

## Formula

\[
\text{Accuracy} =
\frac{\text{Number of Correct Predictions}}
{\text{Total Number of Predictions}}
\]

## Implementation

```python
import numpy as np

def accuracy_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Calculate the accuracy score.

    Args:
        y_true: True labels.
        y_pred: Predicted labels.

    Returns:
        Accuracy score as a float.
    """
    correct_predictions = np.sum(y_true == y_pred)
    total_predictions = len(y_true)

    accuracy = correct_predictions / total_predictions

    return accuracy
```

## Example

```python
y_true = np.array([1, 0, 1, 1, 0, 1])
y_pred = np.array([1, 0, 0, 1, 0, 1])

output = accuracy_score(y_true, y_pred)

print(output)
```

## Output

```text
0.8333333333333334
```

## Explanation

The comparison

```python
y_true == y_pred
```

returns:

```text
[ True  True  False  True  True  True ]
```

NumPy treats `True` as `1` and `False` as `0`, so:

```python
np.sum(y_true == y_pred)
```

returns:

```text
5
```

There are `6` predictions in total:

```text
Accuracy = 5 / 6
         = 0.8333333333333334
```

Therefore, the model accuracy is approximately:

```text
83.33%
```

## Requirements

- Python 3
- NumPy

Install NumPy with:

```bash
pip install numpy
```
