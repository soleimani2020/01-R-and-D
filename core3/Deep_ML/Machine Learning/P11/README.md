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


## Example

```python
binomial_probability(5, 2, 0.5)
# 0.3125
```
Video Tutorial: https://www.youtube.com/watch?v=B-RnkT_fbXI
