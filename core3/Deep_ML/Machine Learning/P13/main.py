import numpy as np

def descriptive_statistics(data: list | np.ndarray) -> dict:
    data = np.array(data)

    # Mean, median
    mean = np.mean(data)
    median = np.median(data)

    # Mode
    values, counts = np.unique(data, return_counts=True)
    mode = values[np.argmax(counts)]

    # Population variance & standard deviation
    variance = np.var(data)      # population
    std_dev = np.std(data)       # population

    # Percentiles
    p25 = np.percentile(data, 25)
    p50 = np.percentile(data, 50)
    p75 = np.percentile(data, 75)
    iqr = p75 - p25

    return {
        "mean": round(mean, 4),
        "median": round(median, 4),
        "mode": mode,
        "variance": round(variance, 4),
        "standard_deviation": round(std_dev, 4),
        "25th_percentile": round(p25, 4),
        "50th_percentile": round(p50, 4),
        "75th_percentile": round(p75, 4),
        "interquartile_range": round(iqr, 4),
    }
