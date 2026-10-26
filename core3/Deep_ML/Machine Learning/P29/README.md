# Linear Regression Using the Normal Equation

> 📊 **Machine Learning** · 📈 **Linear Regression** · 🧮 **Linear Algebra** · 🎯 **Least-Squares Optimization**

## 📌 Project Overview

The **Normal Equation** is an analytical method for finding the coefficients of a linear regression model by minimizing the sum of squared prediction errors.

Unlike gradient descent, which iteratively updates the model parameters, the Normal Equation provides a direct solution when the required matrix operations are well-defined.

This guide derives the Normal Equation step by step using matrix algebra and differentiation.

### 🎥 Video Tutorial

[Linear Regression Using the Normal Equation](https://www.youtube.com/watch?v=1kkVEcmhkL8)

## 🎯 Main Goal

Understand how to derive the least-squares solution for linear regression by minimizing the squared error between the true target values and the model's predictions.

## 🧮 Mathematical Derivation

### 1. Define the Linear Regression Model

We assume the model:

$$
\hat{y} = X\theta
$$

where:

- \(X\) is the feature matrix.
- \(\theta\) is the vector of model coefficients.
- \(y\) is the vector of true target values.
- \(\hat{y}\) is the vector of predicted target values.

If the model includes an intercept, a column of ones must be added to \(X\).

### 2. Define the Residuals and Cost Function

The residual vector represents the difference between the actual and predicted values:

$$
e = y - X\theta
$$

We minimize the sum of squared errors:

$$
J(\theta) = \|y-X\theta\|^2
$$

Using the squared-norm identity:

$$
J(\theta) = (y-X\theta)^\top(y-X\theta)
$$

### 3. Expand the Cost Function

Expanding the matrix product gives:

$$
J(\theta)
= y^\top y
- y^\top X\theta
- \theta^\top X^\top y
+ \theta^\top X^\top X\theta
$$

The terms \(y^\top X\theta\) and \(\theta^\top X^\top y\) are scalars. A scalar is equal to its transpose, so:

$$
y^\top X\theta = \theta^\top X^\top y
$$

Therefore, the cost function simplifies to:

$$
J(\theta)
= y^\top y
- 2\theta^\top X^\top y
+ \theta^\top X^\top X\theta
$$

### 4. Differentiate with Respect to the Parameters

We differentiate each term with respect to \(\theta\).

The first term is independent of \(\theta\):

$$
\frac{\partial}{\partial\theta}(y^\top y)=0
$$

The derivative of the second term is:

$$
\frac{\partial}{\partial\theta}
(-2\theta^\top X^\top y)
= -2X^\top y
$$

The derivative of the quadratic term is:

$$
\frac{\partial}{\partial\theta}
(\theta^\top X^\top X\theta)
= 2X^\top X\theta
$$

This follows because \(X^\top X\) is symmetric.

Combining the terms gives the gradient:

$$
\nabla_\theta J
= -2X^\top y + 2X^\top X\theta
$$

### 5. Set the Gradient to Zero

At the minimum of the least-squares cost function, the gradient must be zero:

$$
\nabla_\theta J=0
$$

Substituting the gradient:

$$
-2X^\top y+2X^\top X\theta=0
$$

Dividing both sides by 2:

$$
-X^\top y+X^\top X\theta=0
$$

Rearranging gives:

$$
X^\top X\theta=X^\top y
$$

This is the **Normal Equation**.

### 6. Solve for the Coefficients

If \(X^\top X\) is invertible, multiply both sides from the left by its inverse:

$$
(X^\top X)^{-1}X^\top X\theta
=
(X^\top X)^{-1}X^\top y
$$

Since:

$$
(X^\top X)^{-1}X^\top X=I
$$

we obtain:

$$
I\theta=(X^\top X)^{-1}X^\top y
$$

Because \(I\theta=\theta\), the final solution is:

$$
\boxed{\theta=(X^\top X)^{-1}X^\top y}
$$

This is the inverse-based form of the least-squares solution.

## 🐍 Python Implementation

The following example uses NumPy to calculate the coefficients with the Normal Equation.

```python
import numpy as np

# Example feature values
X = np.array([
    [1, 1],
    [1, 2],
    [1, 3],
    [1, 4],
    [1, 5]
], dtype=float)

# Target values
y = np.array([2, 4, 5, 4, 6], dtype=float)

# Solve the least-squares problem
theta = np.linalg.lstsq(X, y, rcond=None)[0]

# Predictions
y_pred = X @ theta

print("Coefficients:", theta)
print("Predictions:", y_pred)
```

The first column of ones represents the intercept. The second column contains the feature values.

`np.linalg.lstsq()` computes a least-squares solution without explicitly calculating the matrix inverse. This is generally preferable numerically to implementing the inverse-based formula directly.

## ⚠️ Important Considerations

### When Does the Inverse-Based Formula Work?

The expression:

$$
\theta=(X^\top X)^{-1}X^\top y
$$

requires \(X^\top X\) to be invertible. This holds when the columns of \(X\) are linearly independent.

If features are perfectly collinear or there are insufficient independent observations, the matrix may be singular.

In such cases, a least-squares solver or the Moore–Penrose pseudoinverse can be used.

### Normal Equation vs. Gradient Descent

| Feature | Normal Equation | Gradient Descent |
|---|---|---|
| Approach | Direct least-squares solution | Iterative optimization |
| Learning rate | Not required | Required |
| Iterations | Not required | Usually required |
| Large feature sets | Matrix operations can be expensive | Can be more practical |
| Convergence tuning | No learning-rate tuning | Learning rate and stopping criteria matter |

The Normal Equation is particularly useful for understanding the analytical solution to linear regression, while iterative optimization methods can be more practical for large-scale problems.

## 🚀 Key Takeaways

- Linear regression predicts targets using a linear combination of features.
- Least squares minimizes the sum of squared residuals.
- Differentiating the cost function and setting its gradient to zero produces the Normal Equation.
- The inverse-based solution requires \(X^\top X\) to be invertible.
- Numerical least-squares solvers are generally preferable to explicitly calculating a matrix inverse.

## 📌 Final Formula

The least-squares coefficients are obtained by solving:

$$
X^\top X\theta=X^\top y
$$

When \(X^\top X\) is invertible, the solution is:

$$
\boxed{\theta=(X^\top X)^{-1}X^\top y}
$$

**Project focus:** Linear model → Squared-error cost function → Matrix expansion → Gradient calculation → Normal Equation → Least-squares solution
