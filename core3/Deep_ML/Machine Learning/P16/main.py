import numpy as np

def gradient_descent(X, y, weights, learning_rate, n_epochs, batch_size=1, method='batch'):
    """
    Perform gradient descent optimization.
    
    Args:
        X: Feature matrix of shape (m, n)
        y: Target values of shape (m,)
        weights: Initial weights of shape (n,)
        learning_rate: Step size for gradient descent
        n_epochs: Number of complete passes through the dataset
        batch_size: Size of batches for mini-batch gradient descent (default: 1)
        method: Type of gradient descent ('batch', 'stochastic', or 'mini_batch')
    
    Returns:
        Optimized weights
    """
    # Your code here
    m , n = X.shape 
    y = y.reshape(-1,1) # (m,1)
    
    if weights is None:
        weights = np.zeros((n,1))
    else:
        weights = weights.reshape(-1,1)
    
    
    for epoch in range(n_epochs):
        if method == "batch":
            yhat = X @ weights # (m,1) = (m,n) | (n,1)
            error = yhat - y   # (m,1)
            grad = (2/m) * (X.T @ error)   # (n,1) = (n,m) | ( m,1)
            weights = weights - learning_rate * grad 
            
        elif method == "stochastic":
            for i in range(m):
                x_i = X[i].reshape(1,-1)
                y_i = y[i]
                error = x_i @ weights - y_i
                grad = 2 * x_i.T * error
                weights = weights - learning_rate * grad 
                
                
        elif method == "mini_batch":
            for i in range(0 , m , batch_size):
                X_batch = X[i : i+batch_size]  #(b,1)
                y_batch = y[i : i+batch_size]  #(b,1)
                b = X_batch.shape[0]   # number of samples in a batch 
                error = X_batch @ weights - y_batch
                grad = (2/b) * (X_batch.T @ error)
                weights = weights - learning_rate * grad
                
    return weights.flatten()
                
                
    
        
        
        
        
        
        
        
        
        
        
        
        
        
