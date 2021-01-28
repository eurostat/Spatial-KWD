# -*- coding: utf-8 -*-
#
# @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
#                via Ferrata, 1, I-27100, Pavia, Italy
#
#  @author stefano.gualandi@gmail.com (Stefano Gualandi)
#

from distutils.core import setup, Extension

from Cython.Build import cythonize

extensions = Extension(
    "KWD", ["histogram2D.pyx"],
    extra_compile_args=[
        '-Wno-unused-function', '-std=c++11', '-fopenmp', '-O2', '-ffast-math',
        '-march=native', '-DNDEBUG', '-fno-wrapv'
    ],
    extra_link_args=['-fopenmp', '-O2', '-lm', '-pthread', '-fno-wrapv'])

with open('README.md', encoding="utf-8") as f:
    long_descr = f.read()

setup(name='Spatial-KWD',
      version='0.2.0',
      description='Spatial KWD for Large Spatial Maps',
      author='Stefano Gualandi',
      author_email='stefano.gualandi@gmail.com',
      url='https://github.com/eurostat/Spatial-KWD',
      platforms=['linux', 'macosx', 'windows'],
      download_url=
      'https://github.com/eurostat/Spatial-KWD/archive/v0.2.0-alpha.tar.gz',
      setup_requires=['numpy>=1.16', 'cython>=0.23'],
      install_requires=['numpy>=1.16'],
      long_description=long_descr,
      long_description_content_type='text/markdown',
      ext_modules=cythonize(extensions,
                            compiler_directives={'language_level': "3"}))
