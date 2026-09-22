# Implement Adam Optimization Algorithm

## Overview

This project implements the **Adam (Adaptive Moment Estimation)** optimization algorithm from scratch in Python.

Adam is a gradient-based optimization algorithm widely used in **machine learning** and **deep learning**. It combines ideas from **Momentum** and **RMSProp** by maintaining estimates of both the first and second moments of the gradients.

The optimizer automatically adapts the effective learning rate for each parameter.

---

## Adam Algorithm

For each iteration \(t\), Adam performs the following steps.

### 1. Compute the gradient

$$
g_t = \nabla f(x_t)
$$

### 2. Update the first moment estimate

The first moment acts similarly to momentum:

$$
m_t = \beta_1 m_{t-1} + (1-\beta_1)g_t
$$

### 3. Update the second moment estimate

The second moment keeps track of the squared gradients:

$$
v_t = \beta_2 v_{t-1} + (1-\beta_2)g_t^2
$$

### 4. Bias correction

Because both moment estimates start at zero, they are biased toward zero during the first iterations.

Adam corrects this bias using:

$$
\hat{m}_t =
\frac{m_t}{1-\beta_1^t}
$$

$$
\hat{v}_t =
\frac{v_t}{1-\beta_2^t}
$$

### 5. Update the parameters

Finally, the parameters are updated as:

$$

x_t =
x_{t-1}
-
\alpha
\frac{\hat{m}_t}
{\sqrt{\hat{v}_t}+\epsilon}

$$

where \(\alpha\) is the learning rate.

---

## Function

The main function is:

```python
adam_optimizer(
    f,
    grad,
    x0,
    learning_rate=0.001,
    beta1=0.9,
    beta2=0.999,
    epsilon=1e-8,
    num_iterations=10
)
```

### Parameters

| Parameter        | Description                            |  Default |
| ---------------- | -------------------------------------- | -------: |
| `f`              | Objective function to minimize         | Required |
| `grad`           | Function computing the gradient of `f` | Required |
| `x0`             | Initial parameter values               | Required |
| `learning_rate`  | Optimization step size                 |  `0.001` |
| `beta1`          | Decay rate for the first moment        |    `0.9` |
| `beta2`          | Decay rate for the second moment       |  `0.999` |
| `epsilon`        | Small value for numerical stability    |   `1e-8` |
| `num_iterations` | Number of optimization iterations      |     `10` |

### Returns

The function returns the optimized parameter values after the specified number of iterations.

---

## Implementation

```python
import numpy as np


def adam_optimizer(
    f,
    grad,
    x0,
    learning_rate=0.001,
    beta1=0.9,
    beta2=0.999,
    epsilon=1e-8,
    num_iterations=10
):
    x = np.array(x0, dtype=float)

    # First moment estimate
    m = np.zeros_like(x)

    # Second moment estimate
    v = np.zeros_like(x)

    for t in range(1, num_iterations + 1):

        # Compute gradient
        g = grad(x)

        # Update biased first moment estimate
        m = beta1 * m + (1 - beta1) * g

        # Update biased second moment estimate
        v = beta2 * v + (1 - beta2) * (g ** 2)

        # Bias correction
        m_hat = m / (1 - beta1 ** t)
        v_hat = v / (1 - beta2 ** t)

        # Parameter update
        x -= learning_rate * m_hat / (np.sqrt(v_hat) + epsilon)

    return x
```

---

## Example

Consider the simple quadratic function:

$$
f(x)=x^2
$$

Its gradient is:

$$
\nabla f(x)=2x
$$

We can minimize it using Adam:

```python
import numpy as np


def f(x):
    return np.sum(x ** 2)


def grad(x):
    return 2 * x


x0 = np.array([5.0])

result = adam_optimizer(
    f,
    grad,
    x0,
    learning_rate=0.1,
    num_iterations=100
)

print("Optimized parameters:", result)
print("Objective value:", f(result))
```

The optimizer should gradually move the parameter toward:

```text
x = 0
```

which is the minimum of \(f(x)=x^2\).

---

## Why Adam?

Adam is popular because it combines several useful properties:

* **Adaptive learning rates** for individual parameters
* **Momentum** through the first-moment estimate
* **Gradient scaling** through the second-moment estimate
* **Bias correction** during early iterations
* Good performance for large and noisy optimization problems

Adam is especially common when training neural networks.

---

## Default Hyperparameters

The commonly used Adam hyperparameters are:

```python
learning_rate = 0.001
beta1 = 0.9
beta2 = 0.999
epsilon = 1e-8
```

These defaults work reasonably well for many machine-learning problems, although they may need tuning depending on the objective function.

---

## Requirements

* Python 3
* NumPy

Install NumPy with:

```bash
pip install numpy
```

---

## Key Concepts

This implementation demonstrates:

* Gradient-based optimization
* Exponential moving averages
* Momentum
* Adaptive learning rates
* Bias correction
* Numerical stability
* NumPy vector operations

---

## License

This project is intended for educational purposes and can be freely modified or extended.
