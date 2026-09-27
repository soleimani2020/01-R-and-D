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

# Linear Regression Using the Normal Equation

We start with

$$
\hat{y} = X\theta
$$

and minimize the squared error:

$$
J(\theta) = \|y - X\theta\|^2
$$

Expand:

$$
J(\theta)
=
(y-X\theta)^T(y-X\theta)
$$

$$
J(\theta)
=
y^Ty
-
2\theta^TX^Ty
+
\theta^TX^TX\theta
$$

Differentiate with respect to $\theta$:

$$
\nabla_\theta J
=
-2X^Ty
+
2X^TX\theta
$$

At the minimum:

$$
\nabla_\theta J = 0
$$

Therefore:

$$
-2X^Ty + 2X^TX\theta = 0
$$

$$
X^TX\theta = X^Ty
$$

Multiply by $(X^TX)^{-1}$:

$$
(X^TX)^{-1}X^TX\theta
=
(X^TX)^{-1}X^Ty
$$

Hence:

$$
\boxed{
\theta = (X^TX)^{-1}X^Ty
}
$$
