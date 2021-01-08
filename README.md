# Algorithms for phase unwrapping in Julia

| **License**                     | **Build Status**                                                | **Code Coverage**                                                   |
|:--------------------------------|:----------------------------------------------------------------|:--------------------------------------------------------------------|
| [![][license-img]][license-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

This package implements the phase unwrapping algorithm of Dennis C. Ghiglia and
Louis A. Romero ([*"Robust two-dimensional weighted and unweighted phase
unwrapping that uses fast transforms and iterative
methods"*](https://doi.org/10.1364/JOSAA.11.000107), J. Opt. Soc. Am. A,
Vol. 11, pp. 107-117, 1994).

This package provides the following improvements with respect to the original
method:

- Improved precision (by about 3 digits) by accurately compute `cos(x) - 1`.
  As a result, the algorithm can be used directly in single precision.

- The original algorithm yields the least square solution up to an undetermined
  additive constant.  This bias is avoided (modulo the period) by the
  implemented method.

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://emmt.github.io/Unwrapping.jl/stable

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://emmt.github.io/Unwrapping.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.org/emmt/Unwrapping.jl.svg?branch=master
[travis-url]: https://travis-ci.org/emmt/Unwrapping.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/emmt/Unwrapping.jl?branch=master
[appveyor-url]: https://ci.appveyor.com/project/emmt/Unwrapping-jl/branch/master

[coveralls-img]: https://coveralls.io/repos/emmt/Unwrapping.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/emmt/Unwrapping.jl?branch=master

[codecov-img]: http://codecov.io/github/emmt/Unwrapping.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/emmt/Unwrapping.jl?branch=master
