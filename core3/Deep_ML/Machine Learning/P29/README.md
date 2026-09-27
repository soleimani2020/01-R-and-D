# Linear Regression Using the Normal Equation

We assume the linear regression model

$$
\hat{y} = X\theta
$$

where:

- $X$ is the feature matrix
- $\theta$ is the vector of coefficients
- $y$ is the true target
- $\hat{y}$ is the predicted target

---

## 1. Define the error

The residual is

$$
e = y - X\theta
$$

We minimize the sum of squared errors:

$$
J(\theta) = \lVert y - X\theta \rVert^2
$$

Using the definition of the squared norm:

$$
J(\theta) = (y - X\theta)^T (y - X\theta)
$$

---

## 2. Expand the cost function

Expand the product:

$$
J(\theta) = y^Ty - y^TX\theta - \theta^TX^Ty + \theta^TX^TX\theta
$$

Since

$$
y^TX\theta
$$

is a scalar, its transpose is equal to itself:

$$
y^TX\theta = \theta^TX^Ty
$$

Therefore:

$$
J(\theta) = y^Ty - 2\theta^TX^Ty + \theta^TX^TX\theta
$$

---

## 3. Differentiate with respect to $\theta$

For the first term:

$$
\frac{\partial}{\partial \theta}(y^Ty) = 0
$$

because it does not contain $\theta$.

For the second term:

$$
\frac{\partial}{\partial \theta} \left( -2\theta^TX^Ty \right) = -2X^Ty
$$

For the third term:

$$
\frac{\partial}{\partial \theta} \left( \theta^TX^TX\theta \right) = 2X^TX\theta
$$

because $X^TX$ is symmetric.

Therefore:

$$
\nabla_\theta J = -2X^Ty + 2X^TX\theta
$$

---

## 4. Find the minimum

At the minimum of the cost function:

$$
\nabla_\theta J = 0
$$

So:

$$
-2X^Ty + 2X^TX\theta = 0
$$

Divide both sides by $2$:

$$
-X^Ty + X^TX\theta = 0
$$

Rearrange:

$$
X^TX\theta = X^Ty
$$

This is the **normal equation**.

---

## 5. Solve for $\theta$

We have

$$
X^TX\theta = X^Ty
$$

Multiply both sides from the left by

$$
(X^TX)^{-1}
$$

to get

$$
(X^TX)^{-1}(X^TX)\theta = (X^TX)^{-1}X^Ty
$$

Since

$$
(X^TX)^{-1}(X^TX) = I
$$

we obtain

$$
I\theta = (X^TX)^{-1}X^Ty
$$

and because

$$
I\theta = \theta
$$

the final result is

$$
\theta = (X^TX)^{-1}X^Ty
$$

---

## Final Formula

$$
\theta = (X^TX)^{-1}X^Ty
$$

Video Tutoril: https://www.youtube.com/watch?v=1kkVEcmhkL8

This formula directly gives the least-squares coefficients of the linear regression model.
