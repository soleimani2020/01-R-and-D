# Linear Regression Using the Normal Equation

## Formula

The normal equation is:

$$
\theta = (X^T X)^{-1} X^T y
$$

where:

- $X$ = feature matrix
- $y$ = target vector
- $X^T$ = transpose of $X$
- $\theta$ = regression coefficients

# Derivation of the Normal Equation

We start with the linear regression model:

$$
\hat{y} = X\theta
$$

where:

- $X$ is the feature matrix
- $\theta$ is the vector of coefficients
- $\hat{y}$ is the predicted target

We want the prediction to be as close as possible to the true target $y$.

## 1. Define the squared error

The cost function is:

$$
J(\theta) = \|y - X\theta\|^2
$$

This can be written as:

$$
J(\theta)
=
(y - X\theta)^T(y - X\theta)
$$

---

## 2. Expand the expression

$$
J(\theta)
=
y^Ty
-
y^TX\theta
-
\theta^TX^Ty
+
\theta^TX^TX\theta
$$

The two middle terms are equal because they are scalars:

$$
y^TX\theta
=
\theta^TX^Ty
$$

Therefore:

$$
J(\theta)
=
y^Ty
-
2\theta^TX^Ty
+
\theta^TX^TX\theta
$$

---

## 3. Differentiate with respect to $\theta$

We use:

$$
\frac{\partial}{\partial \theta}(y^Ty)=0
$$

$$
\frac{\partial}{\partial \theta}
\left(
-2\theta^TX^Ty
\right)
=
-2X^Ty
$$

and because $X^TX$ is symmetric:

$$
\frac{\partial}{\partial \theta}
\left(
\theta^TX^TX\theta
\right)
=
2X^TX\theta
$$

Therefore:

$$
\nabla_\theta J
=
-2X^Ty
+
2X^TX\theta
$$

---

## 4. Set the gradient equal to zero

At the minimum:

$$
\nabla_\theta J = 0
$$

Therefore:

$$
-2X^Ty
+
2X^TX\theta
=
0
$$

Divide by $2$:

$$
-X^Ty
+
X^TX\theta
=
0
$$

So:

$$
X^TX\theta
=
X^Ty
$$

This is called the **normal equation**.

---

## 5. Solve for $\theta$

Multiply both sides by $(X^TX)^{-1}$:

$$
(X^TX)^{-1}X^TX\theta
=
(X^TX)^{-1}X^Ty
$$

Since:

$$
(X^TX)^{-1}(X^TX)=I
$$

we obtain:

$$
I\theta
=
(X^TX)^{-1}X^Ty
$$

and therefore:

$$
\boxed{
\theta
=
(X^TX)^{-1}X^Ty
}
$$

---

## Final Formula

$$
\boxed{
\theta
=
(X^TX)^{-1}X^Ty
}
$$

This formula gives the linear regression coefficients directly.

If $X^TX$ is not invertible, the pseudoinverse can be used instead:

$$
\theta = X^+y
$$
