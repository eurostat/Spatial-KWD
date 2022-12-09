# Changelog

## [2022-12-09] Changelog from v0.4.1:

- Fix a compilation warning for CRAN (remove Italian letters with accent)


## [2021-05-07] Changelog from v0.4.0:

- New feature: support for focus areas.
- Fix a compilation error under Solaris 10. Ready for resubmission to CRAN
- Misc bug fixes.


## [2021-03-06] Changelog from v0.3.0:

- New feature: solve unbalanced optimal transport problems. See parameters "unbalanced" and "unbalanced_cost".
- New feature: solve optimal transport problem with or without precomputing the convex hull. See parameter "convex".
-Misc bug fixes.


## [2021-02-12] Changelog from v0.2.5:

- Fixing issues to the R wrapper raised by CRAN.
- Changing the "recode" option to TRUE by default (R wrapper)


## [2021-02-08] Changelog from v0.2.1 to v0.2.4:

- Fighting with the python wrapper and PyPI.


## [2021-01-28] Changelog v0.2.0:

- Fix several bugs to get ready for submission of the R wrapper to the CRAN repository


## [2021-01-13] Changelog v0.1.5:

- Improved Python wrapper with support to helper functions: compareOneToOne, compareOneToMany, compareAll
- New python example with numpy.array
- Docstring for python
- Bug fixing


## [2021-01-13] Changelog v0.1.4:

- Added three main helper functions
- Significantly improved the R wrapper, along with the documentation
- Added the column generation algorithm, with support for warm restarts
- Bug fixing


## [2020-01-01] Changelog v0.1.3:

- Internal private release for debugging


## [2020-12-05] Changelog v0.1.2:

- Added standalone CLI for the C++ library
- Fix a severe error in the computation of GCD (due to back-compatibility to C++11 instead of C++17)
- Added exact solver via complete bipartite graph
- Added log of progress via parameter
- Added support for internal runtime (to be compared with overall runtime)


## [2020-11-02] Changelog v0.1.1:

- new R example using "data.frame" objects
- fixed errors in c++ code
- added basic support for OpenMP (code is not yet full parallel)
