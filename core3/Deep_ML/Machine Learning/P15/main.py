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
