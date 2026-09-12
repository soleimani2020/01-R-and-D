# Calculate Covariance Matrix

**Difficulty:** Easy
**Topic:** Statistics

Write a Python function to calculate the **covariance matrix** for a given set of vectors.

The input is a list of lists, where each inner list represents a feature and contains its observations.

## Formula

The covariance between two features \(X\) and \(Y\) is:

$$
\operatorname{Cov}(X,Y)=
\frac{1}{n-1}
\sum_{i=1}^{n}(x_i-\bar{x})(y_i-\bar{y})
$$

where:

* \(n\) = number of observations
* \(\bar{x}\) = mean of feature \(X\)
* \(\bar{y}\) = mean of feature \(Y\)

The function should return the covariance matrix as a **list of lists**.

## Example

**Input:**

```text
[
    [1, 2, 3],
    [2, 4, 6]
]
```

**Output:**

```text
[
    [1.0, 2.0],
    [2.0, 4.0]
]
```

Include test cases to verify that the implementation works correctly.
