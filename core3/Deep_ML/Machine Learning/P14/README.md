# Calculate Covariance Matrix

**Difficulty:** Easy
**Category:** Statistics

## Description

Implement a Python function that calculates the **covariance matrix** for a given dataset.

The input is a **list of lists**, where each inner list represents a feature and contains its observations.

The function should return the resulting covariance matrix as a **list of lists**.

## Covariance Formula

For two features \(X\) and \(Y\), the sample covariance is:

$$
\operatorname{Cov}(X,Y)
=
\frac{1}{n-1}
\sum_{i=1}^{n}
(x_i-\bar{x})(y_i-\bar{y})
$$

where:

* \(n\) is the number of observations
* \(\bar{x}\) is the mean of \(X\)
* \(\bar{y}\) is the mean of \(Y\)

The diagonal elements of the covariance matrix represent the **variance of each feature**.

## Example

**Input**

```text id="3nxx47"
[
    [1, 2, 3],
    [2, 4, 6]
]
```

**Output**

```text id="c7sp4s"
[
    [1.0, 2.0],
    [2.0, 4.0]
]
```

## Requirements

* Calculate the mean of each feature.
* Calculate the covariance between every pair of features.
* Return the covariance matrix as a list of lists.
* Include test cases to verify the implementation.
