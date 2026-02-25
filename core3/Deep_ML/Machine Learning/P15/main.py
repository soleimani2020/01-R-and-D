import numpy as np

def linear_regression_gradient_descent(X: np.ndarray, y: np.ndarray, alpha: float, iterations: int) -> np.ndarray:
    """
    Perform linear regression using gradient descent (manual loop version).

    Args:
        X: Feature matrix of shape (m, n), first column should be ones for intercept
        y: Target vector of shape (m,)
        alpha: Learning rate
        iterations: Number of gradient descent iterations
    
    Returns:
        Learned weights as a 1D array of shape (n,)
    """

    m, n = X.shape                 # m = number of samples, n = number of features
    y = y.reshape(-1, 1)           # make y a column vector (m,1)
    theta = np.zeros((n, 1))       # initialize weights (n,1), one per feature

    # Loop over number of iterations
    for _ in range(iterations):
        dLdtheta = np.zeros((n, 1))  # gradient vector initialized to zero

        # Loop over each sample to accumulate gradient
        for xi, yi in zip(X, y):
            xi = xi.reshape(-1, 1)      # make xi a column vector (n,1)
            prediction = xi.T @ theta    # compute prediction for this sample
            error = prediction - yi      # compute scalar error
            dLdtheta += 2 * error * xi   # add gradient contribution for this sample

        # Update weights using average gradient
        theta = theta - alpha * (1/m) * dLdtheta

    return theta.flatten()  # return as 1D array


# --------------------------
# Example usage
# --------------------------
X = np.array([
    [1, 1],  # 1 = intercept, 1 = feature value
    [1, 2],
    [1, 3]
])
y = np.array([2, 4, 6])

theta = linear_regression_gradient_descent(X, y, alpha=0.1, iterations=100)
print("Learned weights (theta):", theta)

###### Dradient dscent 1D

import numpy as np 

np.random.seed(0)

x = np.random.randn(10,1)
y = 2*x + np.random.randn(10,1)

w = 0.0 
b = 0.0 
learning_rate = 0.01

def descent(x, y, w, b, learning_rate):
    dldw = 0.0
    dldb = 0.0 
    N = x.shape[0]
    
    for xi, yi in zip(x, y):
        yhat_i = w * xi + b
        error = yhat_i - yi
        
        dldw += 2 * xi * error
        dldb += 2 * error
        
    dldw /= N
    dldb /= N
    
    w = w - learning_rate * dldw
    b = b - learning_rate * dldb 
    
    return w, b
        
for epoch in range(400):
    
    w, b = descent(x, y, w, b, learning_rate) 
    
    yhat = w * x + b
    loss = np.mean((y - yhat)**2)
    
    if epoch % 50 == 0:
        print(f"Epoch {epoch}, Loss: {loss}")

print("Estimated w:", w)
print("Estimated b:", b)



