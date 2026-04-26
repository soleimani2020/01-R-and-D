import numpy as np 
import math 

def adaboost_fit(X, y, n_clf):
    n_samples, n_features = np.shape(X)
    # w = np.full(n_samples, (1 / n_samples))
    # w = np.ones(n_samples) / n_samples
    w = np.array([1 / n_samples for _ in range(n_samples)])
    clfs = []
    
    for _ in range(n_clf):
        # Inside that loop, add variables to keep track of the best stump:
        # best_error / best_feature / best_threshold / best_polarity
        best_error = np.inf; 
        best_feature = None
        best_threshold = None
        best_polarity = None
        # each feature (1D array) as a candidate split axis
        for i in range(n_features):
            # Isolate the feature
            X_i = X[:, i]
            # Candidate thresholds
            Threshold = np.unique(X_i) # no duplicates 
            prediction = np.zeros(len(X_i))
            for j in range(len(Threshold)):
                current_threshold = Threshold[j]
            
                for polarity in [1, -1]:
                    prediction = np.ones(len(X_i))
                    mismatch = np.ones(len(prediction))
                    for k in range(len(X_i)):
                        if polarity == 1:
                            if X_i[k] < current_threshold:
                                prediction[k] = -1
                            else:
                                prediction[k] = 1
                        else:
                            if X_i[k] > current_threshold:
                                prediction[k] = -1
                            else:
                                prediction[k] = 1
                                
                                
                    mismatch = (prediction != y)         # mismatch = [False, True, True, False]
                    error = np.sum(w * mismatch)         # Error:True(1)
                    
                    if error < best_error:
                        best_error = error
                        best_feature = i 
                        best_threshold = current_threshold
                        best_polarity = polarity
                        best_prediction = prediction.copy()
                        
        # compute alpha
        eps = 1e-10
        alpha = 0.5 * math.log((1 - best_error) / (best_error + eps))
        
        # store classifier
        clfs.append({
            "feature": best_feature,
            "threshold": best_threshold,
            "polarity": best_polarity,
            "alpha": alpha
        })
        
        # update weights & normalization
        # (1) : Increase weights for misclassified points
        # (2) : Decrease weights for correctly classified points
        
        w *= np.exp(-alpha * y * best_prediction)
        w /= np.sum(w) 
    
    return clfs
        
                        
                        
X = np.array([[1, 2], [2, 3], [3, 4], [4, 5]])
y = np.array([1, 1, -1, -1])
n_clf = 3

clfs = adaboost_fit(X, y, n_clf)
print(clfs)                    
                        

                        
                                
         

            

   
