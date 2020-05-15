# lap2d
A fast direct dense solver with machine accuracy for 2-D Laplace's equation.

## Overview
By classical potential theory, Laplace's equation can be reformulated as Fredholm integral equation of the second kind on the boundary of the domain, which leads to both a reduction of dimensionality and a bounded operator. However, the discretization of the boundary integral equation will result in a dense matrix. Fortunately, the resulting matrix is rank-structured, to be more precise, hierarchical semi-separable, which allows one to implement an O(n) direct solver for the boundary integral equation.

## Features
* **Machine precision**: Unlike the unbounded Laplacian operator that appeared in the PDE, the formulation of the boundary integral equation will produce a well-conditioned system that is independent of the size of discretization.
* **Direct solver**: Useful for problems involving multiple right-hand sides. Allow for rapid and accurate solutions to relatively ill-conditioned problems.
* **Spectral convergence**: One can prove that the rate of convergence is determined by the order of quadrature rule used when solving the integral equation. By default, the composite 12-point Gauss-Legendre quadrature rule is used.
* **Reduction of dimensionality**: Reduce the problem domain from a two-dimensional one to its boundary (one-dimensional).
* **O(n) time complexity**: Besides low asymptotic cost, the constant of the time complexity is also small.
* **Near-boundary potential evaluation**: The function lapfparam_adap allows one to evaluate the near-boundary potential accurately through Legendre expansion and adaptive Gaussian quadrature after the density of charges/dipoles is solved.

## Caveats
* lap2d currently only supports interior Dirichlet problems, although the boundary integral equation will work for both interior/exterior and Dirichlet/Neumann problems.
* The number of panels in the discretization must be a power of 2. (But one is allowed to adjust the number of Gauss-Legendre nodes on each panel.)
* lap2d depends on NumericalToolbox written by [Serkh](http://www.math.toronto.edu/~kserkh/) et al. The link to the library will be updated here after publication.

## Acknowledgements
I would like to thank my advisor, Professor [Kirill Serkh](http://www.math.toronto.edu/~kserkh/), for both his huge amount of effort devoted and his seemingly infinite amount of patience, support and guidance. I'm also grateful for the opportunity to work on this fun project.

## Reference
P.G. Martinsson, [Boundary Integral Equations and the Nyström method](https://amath.colorado.edu/faculty/martinss/2014_CBMS/Lectures/lecture08.pdf), CBMS/NSF Conference on Fast Direct Solvers (2014).

A. Gillman, P. Young, and P.G. Martinsson, [A direct solver with O(N) complexity for integral equations on one-dimensional domains](https://arxiv.org/pdf/1105.5372.pdf), Front. Math. China, 7 (2012), pp. 217–247.

J. Levandosky, [Notes on potential theory](https://web.stanford.edu/class/math220b/handouts/potential.pdf) (2003).

