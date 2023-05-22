# Example documentation

## Optimize on Circle
In this example, we are interested in find a point on the unit disc that minimizes a two-dimensional convex quadratic function.

If the objective function is not symmetric and the unconstrained minimizer lies within the interior of the unit circle, then the problem has two local minima. With this motivation, the problem defines a nice toy problem in order to test the behavior of the algorithm.

The original problem does not directly fit into the LCQPow framework. In order to fit describe the problem in the LCQPow language, we discretize the unit disc with $N$ points each of which is described by 

```math
x_i = \begin{pmatrix}
    \cos(2 \pi i/N) \\
    \sin(2\pi i/N)
\end{pmatrix}
```

for $i = 1,2,\dots,N$. The convex hull of these points defines a set that approximates the closed unit disc for $N \to \infty$. Let us introduce the slack variable $\lambda_i$ for the constraint 

```math
\begin{pmatrix}
    \cos(2 \pi i/N) \\
    \sin(2\pi i/N)
\end{pmatrix}^\top x + \lambda_i = 1.
```

Then $\lambda_i=0$ implies that $x = x_i$. 

Finally, we introduce a convex combination $ \sum_{i=1}^N \theta_i = 1,$ where each $\theta_i \geq 0$. Now notice that if the vectors of $\theta$ and $\lambda$ satisfy the complementarity constraint

```math
0 \leq \lambda \perp \theta \geq 0,
```
i.e., if each $\lambda_i \theta_i = 0$, then exactly one slack variable $\lambda_{j} = 0$ and $x = x_j$ for some $j$ where $\theta_j = 1$. If complementarity is violated, then the point $x$ lies in the interior of the unit disc. The algorithm thus *converges* to the boundary of the unit disc from its interior for growing complementarity satisfaction.
