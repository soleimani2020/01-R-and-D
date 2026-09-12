# Calculate Covariance Matrix

**Difficulty:** Easy
**Category:** Statistics

## Description

Implement a Python function that calculates the **covariance matrix** for a given dataset.

The input is a list of lists, where each inner list represents a **feature** and contains its observations.

The function should return the covariance matrix as a **list of lists**.

## Covariance Formula

The sample covariance between two features `X` and `Y` is:

$$
\mathrm{Cov}(X,Y) = \frac{1}{n-1}\sum_{i=1}^{n}(x_i-\bar{x})(y_i-\bar{y})
$$

where:

* `n` = number of observations
* `x̄` = mean of feature X
* `ȳ` = mean of feature Y
* `xᵢ`, `yᵢ` = individual observations

The diagonal elements of the covariance matrix represent the **variance of each feature**.

## Example

**Input:**

```text id="6rb0xi"
[
    [1, 2, 3],
    [2, 4, 6]
]
```

**Output:**

```text id="ux4vkr"
[
    [1.0, 2.0],
    [2.0, 4.0]
]
```

## Requirements

* Calculate the mean of each feature.
* Calculate the covariance between every pair of features.
* Construct the covariance matrix.
* Return the result as a list of lists.
* Include test cases to verify the implementation.

* Video Tutorial: https://www.youtube.com/watch?v=WBlnwvjfMtQ
