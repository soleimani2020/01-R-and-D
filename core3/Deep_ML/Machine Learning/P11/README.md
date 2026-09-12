# Binomial Distribution Probability

**Difficulty:** Medium
**Topic:** Probability

Calculate the probability of getting exactly `k` successes in `n` independent Bernoulli trials, where each trial has probability `p` of success.

## Formula

$$
P(X=k)=\binom{n}{k}p^k(1-p)^{n-k}
$$

where:

* `n` = number of trials
* `k` = number of successes
* `p` = probability of success

## Python Solution

```python
import math

def binomial_probability(n, k, p):
    val = math.comb(n, k) * p**k * (1 - p)**(n - k)
    return round(val, 5)
```

## Example

```python
binomial_probability(5, 2, 0.5)
# 0.3125
```
