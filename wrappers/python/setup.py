# -*- coding: utf-8 -*-
#
# @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
#                via Ferrata, 1, I-27100, Pavia, Italy
#
#  @author stefano.gualandi@gmail.com (Stefano Gualandi)
#

from setuptools import setup, Extension, Command, find_packages

from Cython.Build import cythonize

with open('README.md', encoding="utf-8") as f:
    long_descr = f.read()

extensions = Extension(
    name="KWD",
    sources=["KWD_NetSimplex.h", "KWD_Histogram2D.h", "histogram2D.pyx"],
    extra_compile_args=[
        '-Wno-unused-function', '-std=c++11', '-fopenmp', '-O2', '-ffast-math',
        '-march=native', '-DNDEBUG', '-fno-wrapv'
    ],
    extra_link_args=['-fopenmp', '-O2', '-lm', '-pthread', '-fno-wrapv'])

setup(name='Spatial-KWD',
      version='0.2.1',
      packages=find_packages(),
      description='Spatial KWD for Large Spatial Maps',
      author='Stefano Gualandi',
      author_email='stefano.gualandi@gmail.com',
      url='https://github.com/eurostat/Spatial-KWD',
      platforms=['linux', 'macosx', 'windows'],
      download_url=
      'https://github.com/eurostat/Spatial-KWD/archive/v0.2.1-alpha.tar.gz',
      setup_requires=['numpy', 'cython'],
      install_requires=['numpy'],
      long_description=long_descr,
      long_description_content_type='text/markdown',
      ext_modules=cythonize(extensions,
                            compiler_directives={'language_level': "3"}))
