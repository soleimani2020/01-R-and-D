# Poisson Distribution Probability Calculator

**Difficulty:** Easy
**Topic:** Probability

Calculate the probability of observing exactly `k` events using the **Poisson distribution**:

$$
P(X=k)=\frac{e^{-\lambda}\lambda^k}{k!}
$$

where:

* `k` = number of events
* `lam` = mean number of events (λ)

The result is rounded to **5 decimal places**.

## Python Solution

```python
import math

def poisson_probability(k: int, lam: float) -> float:
    probability = (math.exp(-lam) * lam**k) / math.factorial(k)
    return round(probability, 5)
```

## Example

```python
poisson_probability(3, 2)
# 0.18045
```
