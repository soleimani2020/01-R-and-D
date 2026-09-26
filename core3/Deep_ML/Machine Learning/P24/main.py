import numpy as np

def r_squared(y_true, y_pred):
    ss_fit = np.sum((y_true - y_pred) ** 2)
    ss_mean = np.sum((y_true - y_true.mean()) ** 2)

    return 1 - (ss_fit / ss_mean)


y_true = np.array([1, 2, 3, 4, 5])
y_pred = np.array([1.1, 2.1, 2.9, 4.2, 4.8])

print(r_squared(y_true, y_pred))
