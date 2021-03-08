# Spatial-KWD

## <a name="Description"></a>Computing Kantorovich-Wasserstein distances for large spatial maps

This software library contains efficient implementations of Discrete Optimal Transport algorithms for the computation of Kantorovich-Wassestein distances [1] customized for large spatial maps.

The core library is written in standard ANSI-C++11, but it has two wrappers:

1. A [Python wrapper](https://github.com/eurostat/Spatial-KWD/tree/main/wrappers/python)
2. An [R wrapper](https://github.com/eurostat/Spatial-KWD/tree/main/wrappers/R)

If you need a wrapper in another language, please let us know by posting a request on GitHub.

Currently, the Spatial-KWD library is tested on

* Google Colabs notebooks
* Windows 10 (R v4.0.3, Python >= 3.6)
* Mac OS X Bug Sur 11.0.1 (R v4.0.3, Python 3.8.3)
* Linux 20.04.1 LTS (R v4.0.0, Python 3.8.5)

## <a name="Requirements"></a>Basic Usage: Colab Notebooks

The simplest way to test this library is to run one of the following notebooks on Colab:

| Data | Notebook | Link |
|:-|:-|:-|
|**[2021/02/20]**|*Tutorial 1: Using Spatial-KWD with Python*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/Spatial-KWD/blob/main/notebooks/Spatial_KWD_Tutorial_1.ipynb)|
|**[2020/11/22]**|*Tutorial 2: Using Spatial-KWD with R*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/Spatial-KWD/blob/main/notebooks/Spatial_KWD_with_R_Tutorial_2.ipynb)|


## <a name="Python-wrapper"></a>Python Wrapper

To compile the [Python](https://www.python.org/) wrapper you need the following library:

* [Cython](https://cython.org/) (>= 0.23)

To install Cython you can look at the official documentation [Installing Cython](https://cython.readthedocs.io/en/latest/src/quickstart/install.html).

Then, you have to download this repository, move to the subdirectory `.\wrapper\python\`, and then run from command line the following command:

```
python setup.py build_ext
```

This will compile the C++ code and build the python wrapper. If you prefer to compile the wrapper in a local directory, you can run the command:

```
python setup.py build_ext --inplace
```

In this case, however, you can only use the library in the local directory where you compile it.

Once you have done with the installation of the python wrapper, you can test it using one of the following script:

* [test_kwd.py](https://github.com/eurostat/Spatial-KWD/blob/main/examples/test_kwd.py): A basic example with three histograms.
* [test_pandas.py](https://github.com/eurostat/Spatial-KWD/blob/main/examples/test_pandas.py): An example based on Pandas Series, which uses also an exact solver.


## <a name="R-wrapper"></a>R Wrapper

The wrapper for R is contained in the [SpatialKWD_0.3.0.tar.gz](https://github.com/eurostat/Spatial-KWD/releases/download/v0.3.0-alpha/SpatialKWD_0.3.0.tar.gz) package, and all the source code are freely available on this site under the directory [wrappers/R](https://github.com/eurostat/Spatial-KWD/tree/main/wrappers/R).

For compiling the wrapper, you need a recent version of [R](https://www.r-project.org/) and the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) package (we test it with Rcpp v.1.0.5, but it should work also with older versions).

Windows users can download a pre-compiled binary package at the following link:

* [SpatialKWD_0.3.0.zip](https://github.com/eurostat/Spatial-KWD/releases/download/v0.3.0-alpha/SpatialKWD_0.3.0.zip)


Once you have installed the Spatial-KWD package, you can test it running one of the following scripts:

* [test_SKWD.R](https://github.com/eurostat/Spatial-KWD/blob/main/examples/test_SKWD.R): to compute the distances among three dummies histograms
* [test_SKWD_dataframe.R](https://github.com/eurostat/Spatial-KWD/blob/main/examples/test_SKWD_dataframe.R): to compute the distance between two histograms stored in a `data.frame` or in a `data.table`.

If you need an interface for a different R data structure, again, please let us know.


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
