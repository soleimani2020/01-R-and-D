# Linear Regression Using the Normal Equation

> 📈 A mathematical approach to linear regression using **linear algebra, matrix operations, and analytical optimization**.

## 📌 Overview

Linear regression is a supervised machine learning algorithm used to model the relationship between input features and a continuous target variable.

The **Normal Equation** provides an analytical solution for estimating the model parameters by minimizing the sum of squared prediction errors. Unlike gradient descent, it does not require iterative parameter updates or a learning rate.

📚 **Reference:** [Linear Regression Using the Normal Equation](https://www.youtube.com/watch?v=9sK7p2f6Y0E)

## 🎯 Objectives

- Understand the mathematical foundations of linear regression.
- Express the model using matrix notation.
- Derive the cost function and its gradient.
- Obtain the Normal Equation by setting the gradient to zero.
- Understand the advantages and limitations of the analytical solution.
- Compare the Normal Equation with gradient descent.

## 🧠 Mathematical Background

### 1. Linear Regression Model

For a dataset containing \(m\) observations and \(n\) features, the prediction for observation \(i\) is:

$$
\hat{y}^{(i)} = \theta_0 + \theta_1 x_1^{(i)} + \cdots + \theta_n x_n^{(i)}
$$

Where:

- \(\hat{y}^{(i)}\) is the predicted value.
- \(x_j^{(i)}\) is feature \(j\) for observation \(i\).
- \(\theta_0\) is the intercept.
- \(\theta_j\) represents the coefficient associated with feature \(j\).

In matrix notation, the model becomes:

$$
\hat{y} = X\theta
$$

Here, \(X\) is the feature matrix and \(\theta\) is the parameter vector. To include the intercept, the feature matrix contains a column of ones.

### 2. Residuals and Cost Function

A residual is the difference between the actual target value and the predicted value:

$$
e = y - X\theta
$$

The least-squares cost function minimizes the sum of squared residuals:

$$
J(\theta) = \|y-X\theta\|_2^2
$$

Equivalently:

$$
J(\theta) = (y-X\theta)^T(y-X\theta)
$$

Expanding the expression gives:

$$
J(\theta) = y^Ty - 2\theta^TX^Ty + \theta^TX^TX\theta
$$

The objective is to find the parameter vector \(\theta\) that minimizes this cost function.

### 3. Deriving the Normal Equation

To minimize the cost function, we calculate its gradient with respect to the parameter vector:

$$
\nabla_\theta J(\theta)
=
-2X^Ty + 2X^TX\theta
$$

At the minimum, the gradient is zero:

$$
-2X^Ty + 2X^TX\theta = 0
$$

Rearranging gives the **Normal Equation**:

$$
X^TX\theta = X^Ty
$$

If \(X^TX\) is invertible, the parameters can be calculated directly:

$$
\boxed{\theta = (X^TX)^{-1}X^Ty}
$$

This equation provides the least-squares solution without iterative optimization.

## ⚙️ Important Considerations

### When Is the Solution Unique?

The inverse-based solution requires \(X^TX\) to be invertible. This generally requires the columns of \(X\) to be linearly independent.

If features are redundant or highly dependent, the matrix may be singular or poorly conditioned. In such cases, numerical methods based on least squares or the pseudoinverse can be more appropriate than explicitly calculating a matrix inverse.

### Advantages

- Provides a direct analytical solution under the required conditions.
- Does not require a learning rate.
- Does not require iterative gradient updates.
- Is useful for understanding the mathematical foundations of linear regression.

### Limitations

- Matrix operations can become computationally expensive as the number of features increases.
- The inverse-based formula is not suitable when the required inverse does not exist.
- Numerical instability may arise when the feature matrix is poorly conditioned.
- The method remains sensitive to outliers because it minimizes squared errors.

## 🔄 Normal Equation vs. Gradient Descent

| Feature | Normal Equation | Gradient Descent |
|---|---|---|
| Optimization | Analytical solution | Iterative optimization |
| Learning rate | Not required | Required |
| Iterations | Not required | Required |
| Large feature sets | Can be computationally expensive | May be more practical |
| Matrix inversion | Used in the direct formula | Not required |
| Convergence | Direct least-squares solution under suitable conditions | Depends on optimization settings |

## ✅ Key Takeaways

- Linear regression can be expressed compactly using matrix notation.
- The least-squares objective minimizes the sum of squared residuals.
- Setting the gradient of the cost function to zero leads to the Normal Equation.
- The inverse-based formula provides a direct solution when \(X^TX\) is invertible.
- For practical numerical computation, least-squares solvers are generally preferable to explicitly calculating the matrix inverse.

## 🔬 Scientific Motivation

The Normal Equation demonstrates how a machine learning model can be derived using mathematical optimization and linear algebra. Understanding this derivation helps build a foundation for more advanced methods, including gradient-based optimization, regularized regression, and numerical machine learning.
