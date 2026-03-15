def gradient_descent(X, y, weights=None, learning_rate=0.01, n_epochs=100, batch_size=1, method='batch'):
    """
    Loop-based gradient descent for linear regression.

    Args:
        X: Feature matrix of shape (m, n) (first column can be 1 for intercept)
        y: Target vector of shape (m,)
        weights: Initial weights, shape (n,)
        learning_rate: Step size
        n_epochs: Number of passes over the dataset
        batch_size: Batch size for mini-batch
        method: 'batch', 'stochastic', 'mini_batch'

    Returns:
        Optimized weights (1D array)
    """

    m, n = X.shape
    y = y.reshape(-1, 1)  # ensure column vector

    if weights is None:
        weights = np.zeros((n, 1))
    else:
        weights = weights.reshape(-1, 1)

    for epoch in range(n_epochs):

        if method == 'batch':
            # Loop-based accumulation of gradients over all samples
            grad_accum = np.zeros((n, 1))
            for i in range(m):
                xi = X[i].reshape(-1, 1)
                yi = y[i]
                error = (xi.T @ weights) - yi
                grad_accum += 2 * error * xi
            weights -= learning_rate * (1/m) * grad_accum

        elif method == 'stochastic':
            # Update weights for each sample individually
            for i in range(m):
                xi = X[i].reshape(-1, 1)
                yi = y[i]
                error = (xi.T @ weights) - yi
                grad = 2 * error * xi
                weights -= learning_rate * grad  # no division by m

        elif method == 'mini_batch':
            # Loop over batches
            for i in range(0, m, batch_size):
                X_batch = X[i:i+batch_size]
                y_batch = y[i:i+batch_size]
                batch_m = X_batch.shape[0]

                grad_accum = np.zeros((n, 1))
                for j in range(batch_m):
                    xi = X_batch[j].reshape(-1, 1)
                    yi = y_batch[j]
                    error = (xi.T @ weights) - yi
                    grad_accum += 2 * error * xi
                weights -= learning_rate * (1/batch_m) * grad_accum

        else:
            raise ValueError("method must be 'batch', 'stochastic', or 'mini_batch'")

    return weights.flatten()
