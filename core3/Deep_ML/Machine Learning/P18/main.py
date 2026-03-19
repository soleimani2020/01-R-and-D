import numpy as np

def soft_threshold(w: np.ndarray, threshold: float) -> np.ndarray:
    return np.sign(w) * np.maximum(np.abs(w) - threshold, 0)


def l1_regularization_gradient_descent(
    X: np.ndarray,
    y: np.ndarray,
    alpha: float = 0.1,
    learning_rate: float = 0.01,
    max_iter: int = 1000,
    tol: float = 1e-4
) -> tuple:

    n_samples, n_features = X.shape   
    weights = np.zeros(n_features)
    bias = 0.0
    y = y.ravel()

    for iteration in range(max_iter):
        yhat = X @ weights + bias
        error = yhat - y

        gradient_mse = (2 / n_samples) * (X.T @ error)

        prev_weights = weights.copy()
        w_temp = weights - learning_rate * gradient_mse

        weights = soft_threshold(w_temp, learning_rate * alpha)

        dldb = (2 / n_samples) * np.sum(error)
        bias -= learning_rate * dldb

        if np.linalg.norm(weights - prev_weights) < tol:
            break

    return weights, bias


# ✅ Main function to test
def main():
    # Simple test: y = x
    X = np.array([[1], [2], [3], [4], [5]])
    y = np.array([1, 2, 3, 4, 5])

    weights, bias = l1_regularization_gradient_descent(
        X, y,
        alpha=0.01,
        learning_rate=0.01,
        max_iter=2000
    )

    y_pred = X @ weights + bias
    mse = np.mean((y - y_pred) ** 2)

    print("Weights:", weights)
    print("Bias:", bias)
    print("Predictions:", y_pred)
    print("MSE:", mse)
    print("MSE < 0.1:", mse < 0.1)


if __name__ == "__main__":
    main()
