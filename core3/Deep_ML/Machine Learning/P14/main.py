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



def calculate_covariance_matrix(vectors: list[list[float]]) -> list[list[float]]:

    n_features = len(vectors)
    n_observations = len(vectors[0])

    means = []

    for feature in vectors:
        mean = sum(feature) / n_observations
        means.append(mean)

    covariance_matrix = []

    for i in range(n_features):
        row = []

        for j in range(n_features):
            cov = 0

            for k in range(n_observations):
                cov += (
                    (vectors[i][k] - means[i]) *
                    (vectors[j][k] - means[j])
                )

            cov = cov / (n_observations - 1)
            row.append(cov)

        covariance_matrix.append(row)

    return covariance_matrix
