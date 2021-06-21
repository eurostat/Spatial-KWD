# Spatial-KWD

GitHub actions: <a class="reference external" href="https://github.com/eurostat/Spatial-KWD/actions"><img alt="Build Status" src="https://github.com/eurostat/Spatial-KWD/workflows/build/badge.svg?branch=master&amp;event=push"></a> 

Downloads of the Python wrapper from PyPI:

<a class="reference external" href="https://badge.fury.io/py/Spatial-KWD"><img alt="PyPI version" src="https://badge.fury.io/py/Spatial-KWD.svg"></a> [![Downloads](https://pepy.tech/badge/spatial-kwd)](https://pepy.tech/project/spatial-kwd) [![Downloads](https://pepy.tech/badge/spatial-kwd/month)](https://pepy.tech/project/spatial-kwd) [![Downloads](https://pepy.tech/badge/spatial-kwd/week)](https://pepy.tech/project/spatial-kwd)

Downloads of the R package from CRAN: 

![](http://cranlogs.r-pkg.org/badges/grand-total/SpatialKWD) ![](http://cranlogs.r-pkg.org/badges/last-week/SpatialKWD)


## <a name="Description"></a>Computing Kantorovich-Wasserstein distances for large spatial maps

This software library contains efficient implementations of Discrete Optimal Transport algorithms for the computation of Kantorovich-Wassestein distances [1] customized for large spatial maps.

The core library is written in standard ANSI-C++11, but it has two wrappers:

1. A [Python wrapper](https://pypi.org/project/Spatial-KWD/) available from PyPI
2. An [R wrapper](https://CRAN.R-project.org/package=SpatialKWD) available from CRAN

If you need a wrapper in another language, please let us know by posting a request on GitHub.

Currently, the Spatial-KWD library is tested on

* Google Colabs notebooks
* Windows 10 (R v4.0.5, Python >= 3.6)
* Mac OS X Bug Sur 11.0.1 (R v4.0.5, Python 3.8.3)
* Linux 20.04.1 LTS (R v4.0.0, Python 3.8.5)
* Solaris 10 (R v4.0.5)

## <a name="Requirements"></a>Basic Usage: Colab Notebooks

The simplest way to test this library is to run one of the following notebooks on Colab, which include an example of installation of the library.

| Data | Notebook | Link |
|:-|:-|:-|
|**[2021/05/10]**|*Tutorial 1: Using Spatial-KWD with Python*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/Spatial-KWD/blob/main/notebooks/Spatial_KWD_Tutorial_1.ipynb)|
|**[2020/05/10]**|*Tutorial 2: Using Spatial-KWD with R*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/Spatial-KWD/blob/main/notebooks/Spatial_KWD_with_R_Tutorial_2.ipynb)|


## <a name="Python-wrapper"></a>Python Wrapper

The simplest way to install the python wrapper is to run the following command:

```
pip install Spatial-KWD
```

For basic usage of the library, please, look at the Colab notebok.

To compile the [Python](https://www.python.org/) wrapper directly from the source code, you need the following library:

* [Cython](https://cython.org/) (>= 0.23)

To install Cython you can look at the official documentation [Installing Cython](https://cython.readthedocs.io/en/latest/src/quickstart/install.html).

Then, you have to download this repository, and run from command line the following command:

```
make buildpython
```

This will compile the C++ code and build the python wrapper. If you prefer to compile the wrapper in a local directory, you can run the command:

```
python setup.py build_ext --inplace
```

In this case, however, you can only use the library in the local directory where you compile it.


## <a name="R-wrapper"></a>R Wrapper

The simplest way to install the R wrapper is to run the following command from the R shell:

```
install.packages("SpatialKWD")
```

For basic usage of the library, please, look at the Colab notebok.

To compile the wrapper directly from the source code, you need a recent version of [R](https://www.r-project.org/) and the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package (we test it with Rcpp v.1.0.5, but it should work also with older versions).

Windows users can download a pre-compiled binary package at the following link:

* [SpatialKWD_0.4.0.zip](https://github.com/eurostat/Spatial-KWD/releases/download/v0.4.0-alpha/SpatialKWD_0.4.0.zip)


If you need an interface for a different R data structure, please, drop us an email.


## <a name="Cpp-CLI"></a>Cpp Command Line Tool

By compiling the Cpp source code, it is possible to use the **SpatialKWD** directly from the command line.
The compilation is ruled by a simple Makefile, and is executed by the command:

```
make laptop
```

This build an executable named `solver` in the `bin` folder, which can be tested with the following command, which use a synthetic data from the whole Belgium region:

```
./bin/solver data/L_1000_100.csv
```

The output shold be as follows:

```
data/L_1000_100.csv
start solver
WARNING: the Xs input coordinates are not consecutives integers.
WARNING: the Ys input coordinates are not consecutives integers.
INFO: Recoding the input coordinates to consecutive integers.
INFO: change <verbosity> to info
INFO: change <opt_tolerance> to 0.000001
Approx => 32140: fobj: 0.689998, time: 2.1860, status: Optimal, iter: 167063, arcs: 58987, nodes: 39302
```

### <a name="About"></a>About

<table align="center">
    <tr> <td align="left"><i>Contributors</i></td>
    <td align="left" valign="middle">
<a href="https://github.com/stegua"><img src="https://github.com/stegua.png" width="40"></a>
</td>  </tr>
    <!-- <tr> <td align="left"><i>version</i></td> <td align="left"> </td> </tr> -->
    <tr> <td align="left"><i>Status</i></td> <td align="left">since 2020</td> </tr>
    <tr> <td align="left"><i>License</i></td> <td align="left"><a href="https://joinup.ec.europa.eu/sites/default/files/custom-page/attachment/2020-03/EUPL-1.2%20EN.txt">EUPL v1.2</a><i></i></td> </tr>
</table>

### <a name="References"></a>Main References

* **[1]** Bassetti F., Gualandi S., Veneroni M. [**On the computation of Kantorovich-Wasserstein distances between 2D-histograms by uncapacitated minimum cost flows**](https://epubs.siam.org/doi/abs/10.1137/19M1261195). SIAM J. Optim., 30(3), 2441–2469, 2020. Preprint on arXiv: [1804.00445](https://arxiv.org/abs/1804.00445).


### Note
Windows package of the R wrapper is built with [Win-Builder](https://win-builder.r-project.org/upload.aspx).
