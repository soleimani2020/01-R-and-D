import numpy as np

def calculate_covariance_matrix(vectors: list[list[float]]) -> list[list[float]]:
    vectors = np.array(vectors, dtype=float)
    
    # Treat rows as variables → transpose
    vectors = vectors.T  # now shape = (n_obs, n_vars)
    
    n = vectors.shape[0]  # number of observations
    means = np.mean(vectors, axis=0)
    centered = vectors - means
    cov_matrix = (centered.T @ centered) / (n-1)
    return cov_matrix.tolist()
