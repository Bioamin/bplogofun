import os
import cython_gsl
from setuptools import setup, Extension
from Cython.Build import cythonize
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = f.read()
#bpexact_ext = Extension('bplogofun.exact', ['src/exact.c'])
nsb_ext = Extension('*', ['src/*.pyx'],
                    include_dirs = [cython_gsl.get_include(),
                                    cython_gsl.get_cython_include_dir()],
                    libraries=cython_gsl.get_libraries(),
                    library_dirs=[cython_gsl.get_library_dir()],
                    )


setup(name = "bplogofun",
      install_requires=['Cython', 'CythonGSL', 'statsmodels', 'numpy', 'scipy', 'pandas', 'patsy'],
      packages = ["bplogofun"],
      package_data={'bplogofun': ['eps/Template.eps']},
      entry_points = {
          "console_scripts": ['bplogofun = bplogofun.bplogofun:main']},
      version = "0.1.1",
      description = "Something Something bplogofun",
      long_description=long_description,
      license='GPLv3',
      url = "www.nowhere.com",
      ext_modules=cythonize(nsb_ext),)
