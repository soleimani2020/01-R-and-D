# 📈 Linear Regression Using Gradient Descent (EASY)

**Difficulty:** Easy  
**Category:** Machine Learning  

This project implements Linear Regression from scratch using Batch Gradient Descent in Python with NumPy.

The goal is to minimize the Mean Squared Error (MSE) loss function and learn the optimal model coefficients (weights).

---

## 📌 Problem Description

Write a Python function that performs linear regression using gradient descent.

The function must:

- Accept:
  - `X`: Feature matrix including a column of ones for the intercept  
  - `y`: Target vector  
  - `alpha`: Learning rate  
  - `iterations`: Number of gradient descent steps  

- Return:
  - The learned coefficient vector `theta` as a NumPy array  

---

## 📐 Mathematical Formulation

We minimize the Mean Squared Error (MSE) loss:

L(θ) = (1 / 2m) * Σ (hθ(xᶦ) − yᶦ)²

Where:

- m = number of training samples  
- θ = parameter vector  
- hθ(x) = Xθ = prediction function  

The factor 1/2 simplifies the gradient expression.

---

## 🔄 Gradient Descent Update Rule

Gradient of the loss function:

∇L(θ) = (1/m) Xᵀ (Xθ − y)

Update step:

θ := θ − α ∇L(θ)

Where:

- α = learning rate  

---

## 📊 Input Specifications

| Variable | Shape | Description |
|----------|--------|------------|
| `X` | `(m, n)` | Feature matrix (includes bias column of ones) |
| `y` | `(m,)` | Target values |
| `theta` | `(n,)` | Model weights |

- m = number of training examples  
- n = number of features (including bias)

---

## 🛠 Implementation

```python
import numpy as np

def linear_regression_gradient_descent(X: np.ndarray,
                                       y: np.ndarray,
                                       alpha: float,
                                       iterations: int) -> np.ndarray:
    """
    Perform Linear Regression using Batch Gradient Descent.

    Parameters:
        X (np.ndarray): Feature matrix of shape (m, n)
        y (np.ndarray): Target vector of shape (m,)
        alpha (float): Learning rate
        iterations (int): Number of gradient descent iterations

    Returns:
        np.ndarray: Learned weight vector of shape (n,)
    """

    m, n = X.shape
    theta = np.zeros(n)  # Initialize weights to zero

    for _ in range(iterations):
        predictions = X @ theta              # (m,)
        errors = predictions - y             # (m,)
        gradient = (1/m) * (X.T @ errors)    # (n,)
        theta = theta - alpha * gradient

    return theta

Video Tutorial: https://www.youtube.com/watch?v=U6Z-UjkJZjA&t=470s

Video Tutorial: https://www.youtube.com/@NicholasRenotte/featured 





