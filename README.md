# lap2d
A fast direct dense solver with machine accuracy for 2-D Laplace's equation.

## Overview
By classical potential theory, Laplace's equation can be reformulated as Fredholm integral equation of the second kind on the boundary of the domain, which leads to both a reduction of dimensionality and a bounded operator. However, the discretization of the boundary integral equation will result in a dense matrix. Fortunately, the resulting matrix is rank-structured, to be more precise, hierarchical semi-separable, which allows one to implement an O(n) direct solver for the boundary integral equation.

## Features
* **Machine precision**: Unlike the unbounded Laplacian operator that appeared in the PDE, the formulation of the boundary integral equation will produce a well-conditioned system that is independent of the size of discretization.
* **Direct solver**: Useful for problems involving multiple right-hand sides. Allow for rapid and accurate solutions to relatively ill-conditioned problems.
* **Spectral convergence**: One can prove that the rate of convergence is determined by the order of quadrature rule used when solving the integral equation. By default, the composite 12-point Gauss-Legendre quadrature rule is used.
* **Reduction of dimensionality**: Reduce the problem domain from a two-dimensional one to its boundary (one-dimensional).
* **Fast**: Besides low asymptotic cost, the constant of the time complexity is also small.
* **Accurate evaluation of near-boundary potential**: The function lapfparam_adap allows one to evaluate near-boundary potential accurately through Legendre expansion and adaptive Gaussian quadrature, given the density of charges/dipoles on the boundary.

## Caveats
* lap2d currently only supports interior Dirichlet problems, although the boundary integral equation will work for both interior/exterior and Dirichlet/Neumann problems.
* The number of panels in the discretization must be a power of 2. (But one is allowed to adjust the number of Gauss-Legendre nodes on each panel.)
* lap2d depends on NumericalToolbox written by [Serkh](http://www.math.toronto.edu/~kserkh/) et al. The link to the library will be updated here after publication.
* lap2d uses a binary tree instead of a quadtree to hierarchically partition the boundary, which causes the theoretical asymptotic time complexity to be O(nlogn). However, it turns out that the coefficient of the nlogn term is tiny and for any reasonable size of discretization, the solver exhibits O(n) asymptotic time complexity.

## Acknowledgements
I would like to thank my advisor, Professor [Kirill Serkh](http://www.math.toronto.edu/~kserkh/), for both his huge amount of effort devoted and his seemingly infinite amount of patience, support and guidance. I'm also grateful for the opportunity to work on this fun project.

## Reference
P.G. Martinsson, [Boundary Integral Equations and the Nyström method](https://amath.colorado.edu/faculty/martinss/2014_CBMS/Lectures/lecture08.pdf), CBMS/NSF Conference on Fast Direct Solvers, 2014.

A. Gillman, P. Young, and P.G. Martinsson, [A direct solver with O(N) complexity for integral equations on one-dimensional domains](https://arxiv.org/pdf/1105.5372.pdf), Front. Math. China, 7 (2012), pp. 217–247.

P.G. Martinsson, [Direct Solvers for Integral Equations](https://amath.colorado.edu/faculty/martinss/2014_CBMS/Lectures/lecture09.pdf), CBMS/NSF Conference on Fast Direct Solvers, 2014.

J. Levandosky, [Notes on potential theory](https://web.stanford.edu/class/math220b/handouts/potential.pdf), Stanford University, 2003.

K. L. Ho and L. Greengard, [A fast direct solver for structured linear systems by recursive skeletonization](https://arxiv.org/pdf/1110.3105.pdf), SIAM J. Sci. Comput., vol. 34, no. 5, pp. 2507–2532, 2012.
