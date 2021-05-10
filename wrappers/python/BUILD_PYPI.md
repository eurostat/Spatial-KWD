# How to build the Python wrapper for submission to PyPY

## Under Windows

Under windows, I have created four different enviroments:

* pyenv36, pyenv37, pyenv38, pyenv39

You can activate the enviroment by the command:

```
> conda activate pyenv37
```

For each enviroment, we have to repeat the following steps:

```
> python3 setup.py sdist bdist_wheel
```

and then

```
> python3 -m twine upload --repository pypi dist/*
```

Beaware that once you upload a version for the first time, later you cannot longer modify the README.md file until the next version.

## Under Linux: *.manylinux2010

The build for linux is a dleicate issue and it is performed under a docker image of an old version of CentOS, using the tool `auditwheel`. Before running the following commands, be sure to have installed on the docker container the library: `numpy, Cython, sdist, auditwheel, twine`.

```
% /opt/python/cp36-cp36m/bin/python3 setup.py sdist bdist_wheel
% auditwheel repair dist/Spatial_KWD-0.4.0-cp36-cp36m-linux_x86_64.whl
% rm dist/Spatial_KWD-0.4.0-cp36-cp36m-linux_x86_64.whl
% mv wheelhouse/Spatial_KWD-0.4.0-cp36-cp36m-manylinux2010_x86_64.whl dist/
```

Repeat for every different python environment: *cp36-cp36m,cp37-cp37m,cp38-cp38,cp39-cp39*.

In the end, you can upload all the library to PyPI with the command:

```
> python3 -m twine upload --repository pypi dist/*
```

## Under MAC OS X

*To be completed.*


### Useful links:

- https://levelup.gitconnected.com/how-to-deploy-a-cython-package-to-pypi-8217a6581f09
- https://malramsay.com/post/perils-of-packaging/
- https://packaging.python.org/tutorials/packaging-projects/
