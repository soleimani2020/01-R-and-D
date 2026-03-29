import math
import numpy as np 

def poisson_probability(k, lam):
    val = (lam**k) * (np.exp(-lam)) / math.factorial(k)
    return round(val, 5)
