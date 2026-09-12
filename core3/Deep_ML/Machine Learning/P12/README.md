# Normal Distribution PDF Calculator

**Difficulty:** Medium
**Topic:** Probability

Calculate the **probability density function (PDF)** of a normal distribution for a given value `x`, mean `mu`, and standard deviation `sigma`.

## Formula

$$
f(x)=\frac{1}{\sigma\sqrt{2\pi}}
e^{-\frac{(x-\mu)^2}{2\sigma^2}}
$$

where:

* `x` = value
* `mu` = mean
* `sigma` = standard deviation


## Example

```python id="j4t7m2"
normal_pdf(0, 0, 1)
# 0.39894
```
