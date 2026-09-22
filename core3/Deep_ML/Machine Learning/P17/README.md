# Implement Adam Optimization Algorithm

## Overview

This project implements the **Adam (Adaptive Moment Estimation)** optimization algorithm from scratch in Python.

Adam is a gradient-based optimization algorithm widely used in **machine learning** and **deep learning**. It combines ideas from **Momentum** and **RMSProp** by maintaining estimates of both the first and second moments of the gradients.

The optimizer automatically adapts the effective learning rate for each parameter.

In Adam, When the gradient is small (close to zero), Adam can take relatively larger steps to move efficiently through flat regions. When the gradient is large, it reduces the effective step size, helping to prevent overshooting the minimum.

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

$$
x_t = x_{t-1} - \alpha \frac{\hat{m}_t}{\sqrt{\hat{v}_t} + \epsilon}
$$

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
