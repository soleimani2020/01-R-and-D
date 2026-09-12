import numpy as np
from statistics import mode

def descriptive_statistics(data):
    
    mean = np.mean(data)
    median = np.median(data)
    mode_value = mode(data)
    
    variance = np.var(data)
    std = np.std(data)
    
    q1 = np.percentile(data, 25)
    q2 = np.percentile(data, 50)
    q3 = np.percentile(data, 75)
    
    iqr = q3 - q1
    
    return {
        "mean": mean,
        "median": median,
        "mode": mode_value,
        "variance": variance,
        "standard_deviation": std,
        "25th_percentile": q1,
        "50th_percentile": q2,
        "75th_percentile": q3,
        "interquartile_range": iqr
    }
