Spatial-KWD
===========

Computing Kantorovich-Wasserstein distances for large spatial maps
---

**<a name="Description"></a>Description**

**<a name="Requirements"></a>Requirements**

To compile the [Python](https://www.python.org/) wrapper you only need the following library:

* Cython (>= 0.23)

While for compiling the [R](https://www.r-project.org/) wrapper, you need:

* Rcpp

The KWD library is tested on:

* Google Colabs
* Windows 10 (python >= 3.6)
* Mac OS X Bug Sur 11.0.1 (python 3.8.3)
* Linux 20.04.1 LTS (python 3.8.5)

**<a name="Quick-launch"></a>Quick launch**

For the python wrapper, you have to download this repository, move to the subdirectory `.\wrapper\python\`, and then run from command line the following command:

```
python setup.py build_ext --inplace
```
This will compile the C++ code and build the python wrapper.

For the R wrapper, we recommend to look at Tutoral 2, the second notebook listed below.


**<a name="Usage"></a>Usage**

An example for using this library is given in the following notebook:

| Data | Notebook | Link |
|:-|:-|:-|
|**[2020/11/13]**|*Tutorial 1: Testing the Spatial-KWD library*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/Spatial-KWD/blob/main/notebooks/Spatial_KWD_Tutorial_1.ipynb)|
|**[2020/11/22]**|*Tutorial 2: Using Spatial-KWD with R*|[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/eurostat/Spatial-KWD/blob/main/notebooks/Spatial_KWD_with_R_Tutorial_2.ipynb)|

**<a name="About"></a>About**

<table align="center">
    <tr> <td align="left"><i>contributors</i></td> 
    <td align="left" valign="middle">
<a href="https://github.com/stegua"><img src="https://github.com/stegua.png" width="40"></a>
</td>  </tr> 
    <!-- <tr> <td align="left"><i>version</i></td> <td align="left"> </td> </tr> -->
    <tr> <td align="left"><i>status</i></td> <td align="left">since 2020</td> </tr> 
    <tr> <td align="left"><i>license</i></td> <td align="left"><!-- <a href="https://joinup.ec.europa.eu/sites/default/files/custom-page/attachment/2020-03/EUPL-1.2%20EN.txt">EUPL</a>--> <i></i></td> </tr> 
</table>

**<a name="References"></a>Main References** 

* Bassetti F., Gualandi S., Veneroni M. (2018): [**On the computation of Kantorovich-Wasserstein distances between 2D-histograms by uncapacitated minimum cost flows**](https://epubs.siam.org/doi/abs/10.1137/19M1261195). SIAM J. Optim., 30(3), 2441–2469, 2020. Preprint on arXiv: [1804.00445](https://arxiv.org/abs/1804.00445).
