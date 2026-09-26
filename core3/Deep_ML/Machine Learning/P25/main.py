import numpy as np

def mae(y_true, y_pred):
    N = y_true.size
    MAE = (1 / N) * np.sum(np.abs(y_true - y_pred))
    return MAE


def mae(y_true, y_pred):
    MAE = np.mean(np.abs(y_true - y_pred))
    return MAE
    
    

y_true = np.array([3, -0.5, 2, 7])
y_pred = np.array([2.5, 0.0, 2, 8])

print(mae(y_true, y_pred))
