import math

def binomial_probability(n: int, k: int, p: float) -> float:
    """
    Calculate the probability of exactly k successes in n Bernoulli trials.
    """
    
    val = (
        math.factorial(n) /
        (math.factorial(k) * math.factorial(n-k))
    ) * p**k * (1-p)**(n-k)

    return val
