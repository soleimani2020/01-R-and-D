import numpy as np

def pegasos_kernel_svm(data: np.ndarray, labels: np.ndarray, kernel='linear', lambda_val=0.01, iterations=100, sigma=1.0) -> tuple:
    """
    Train a kernel SVM using the deterministic Pegasos algorithm.
    
    Args:
        data: Training data of shape (n_samples, n_features)
        labels: Labels of shape (n_samples,) with values in {-1, 1}
        kernel: 'linear' or 'rbf'
        lambda_val: Regularization parameter
        iterations: Number of training iterations
        sigma: RBF kernel bandwidth (only used if kernel='rbf')
    
    Returns:
        Tuple of (alphas, bias) where alphas is a list and bias is a float
    """
    # Your code here
    n_samples , n_features = data.shape 
    b = 0 
    alpha = np.zeros(n_samples)
    
    if kernel == 'linear':
        K = data @ data.T 
        
        
    elif kernel == 'rbf':
        diff = data[:, None, :] - data[None, :, :]
        K = np.exp(-np.sum(diff**2, axis=2) / (2 * sigma**2))
        
        
    for t in range(1, iterations+1):
        j = (t - 1) % n_samples 
        eta_t = 1 /(lambda_val*t)
        
        K_ij = K[:,j]
        f_x = np.sum(alpha * labels * K_ij) + b 
        
        margin = labels[j] * f_x
        
        alpha = alpha * (1 - eta_t * lambda_val)
        
        if margin < 1 :
            alpha[j] = alpha[j] + eta_t
            b = b +  eta_t * labels[j]
            
        #else:
            #print("margin satisfied → no corrective update")
            
            
    return (alpha.tolist() , float(b))  
            
            
            

        
        
        
    
        
