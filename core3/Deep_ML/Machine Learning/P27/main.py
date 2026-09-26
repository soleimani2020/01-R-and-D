import numpy as np
from typing import List, Tuple

def k_fold_cross_validation(n_samples: int, k: int = 5, shuffle: bool = True) -> List[Tuple[List[int], List[int]]]:
    """
    Generate train/test index splits for k-fold cross-validation.
    
    Args:
        n_samples: Total number of samples in the dataset
        k: Number of folds (default 5)
        shuffle: Whether to shuffle indices before splitting (default True)
    
    Returns:
        List of (train_indices, test_indices) tuples
    """
    # Your code here
    indices  = np.arange(0,n_samples)
    
    if shuffle:
        np.random.shuffle(indices)
        
    base_size  =  n_samples // k 
    remainder  =  n_samples % k
    
    """
    n_samples = 10 , k = 3
    base_size = 3  | remainder = 1 
    i = 0 → 0 < 1 → append 4
    i = 1 → 1 < 1 is False → append 3
    i = 2 → 2 < 1 is False → append 3
    
    """
    
    fold_sizes = []
    for i in range(k):
        if i < remainder:
            fold_sizes.append(base_size + 1)  
        else:
            fold_sizes.append(base_size)
            
            
    # indices
    # fold_sizes
    
    folds = []
    start = 0 
    for size in fold_sizes:
        end = start + size
        folds.append(indices[start:end])
        start = end 
        
    
    # folds
    splits = []
    
    for i in range(k):
        test = folds[i]
        train = []
    
        for j in range(k):
            if j != i:
                train.append(folds[j])
    
        train = np.concatenate(train)
    
        splits.append((train.tolist() , test.tolist() )) 
        
        
    return splits
        
                    
                
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
