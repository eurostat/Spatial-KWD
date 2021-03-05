# -*- coding: utf-8 -*-
#
# @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
#                via Ferrata, 1, I-27100, Pavia, Italy
#
#  @author stefano.gualandi@gmail.com (Stefano Gualandi)
#

from setuptools import setup, Extension, Command, find_packages

from Cython.Build import cythonize
import numpy as np

import platform

with open('README.md', encoding="utf-8") as f:
    long_descr = f.read()

CC_ARGS = [
    '-Wno-unused-function', '-std=c++11', '-fopenmp', '-O2', '-ffast-math',
    '-march=native', '-DNDEBUG', '-fno-wrapv'
]

LD_ARGS = ['-lgomp', '-O2', '-lm', '-pthread', '-fno-wrapv']

if platform.system() == 'Windows':
    CC_ARGS = []
    LD_ARGS = []

if platform.system() == 'Darwin':
    CC_ARGS = [
        '-std=c++11', '-stdlib=libc++', '-O2', '-ffast-math', '-DNDEBUG'
    ]
    LD_ARGS = ['-O2', '-lm', '-pthread', '-fno-wrapv']

extensions = Extension(name="KWD",
                       sources=["histogram2D.pyx"],
                       include_dirs=['./'],
                       extra_compile_args=CC_ARGS,
                       extra_link_args=LD_ARGS)

setup(name='Spatial-KWD',
      version='0.2.6',
      packages=find_packages(),
      description='Spatial KWD for Large Spatial Maps',
      author='Stefano Gualandi',
      author_email='stefano.gualandi@gmail.com',
      url='https://github.com/eurostat/Spatial-KWD',
      platforms=['linux', 'macosx', 'windows'],
      download_url=
      'https://github.com/eurostat/Spatial-KWD/archive/v0.2.5-alpha.tar.gz',
      setup_requires=['numpy', 'cython'],
      install_requires=['numpy'],
      include_dirs=np.get_include(),
      long_description=long_descr,
      long_description_content_type='text/markdown',
      ext_modules=cythonize(extensions,
                            compiler_directives={'language_level': "3"}))
