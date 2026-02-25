import numpy as np

def linear_regression_gradient_descent(X: np.ndarray, y: np.ndarray, alpha: float, iterations: int) -> np.ndarray:
    """
    Perform linear regression using gradient descent.

    Args:
        X: Feature matrix of shape (m, n) where first column is all ones (for intercept)
        y: Target vector of shape (m,)
        alpha: Learning rate
        iterations: Number of gradient descent iterations
    
    Returns:
        Learned weights as a 1D array of shape (n,)
    """

    m, n = X.shape
    y = y.reshape(-1, 1)  # Ensure y is a column vector
    ### theta must have the same numebr of rows as columns in X
    theta = np.zeros((n, 1))  # Initialize weights to zeros


    # Prediction = (X.T @ Theta)
    # Error = (X.T @ Theta - yi)
    # Loss = (X.T @ Theta - yi)**2 
    # Loss Derivative (dLdtheta) = 2 *(X.T @ Theta - yi)*X 
    
    for _ in range(iterations):
        dLdtheta = np.zeros((n, 1))  # gradient vector
        for xi , yi in zip(X,y):
            xi = xi.reshape(-1, 1)  # make column vector (n,1)
            prediction = xi.T @ theta
            error = prediction - yi
            dLdtheta += 2 * error * xi  # derivative for this sample
        
        theta = theta - alpha * (1/m) * dLdtheta  # update weights
        
    # Your code here: implement gradient descent

    return theta.flatten()
    
    
# Simple dataset
X = np.array([[1, 1], [1, 2], [1, 3]])  # first column = 1 for intercept
y = np.array([2, 4, 6])                  # target

theta = linear_regression_gradient_descent(X, y, alpha=0.1, iterations=100)
print(theta)



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


