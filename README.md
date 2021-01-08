# Algorithms for phase unwrapping in Julia

This package implements the phase unwrapping algorithm of Dennis C. Ghiglia and
Louis A. Romero (*"Robust two-dimensional weighted and unweighted phase
unwrapping that uses fast transforms and iterative methods"*,
J. Opt. Soc. Am. A, Vol. 11, pp. 107-117, 1994).

This package provides the following improvements with respect to the original
method:

- Improved precision (by about 3 digits) by accurately compute `cos(x) - 1`.
  As a result, the algorithm can be used directly in single precision.
