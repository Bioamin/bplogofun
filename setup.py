import os
from setuptools import setup, Extension
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = f.read()

setup(name = "bplogofun",
      install_requires=['statsmodels', 'numpy', 'scipy', 'pandas', 'patsy'],
      packages = ["bplogofun"],
      package_data={'bplogofun': ['eps/Template.eps']},
      entry_points = {
          "console_scripts": ['bplogofun = bplogofun.bplogofun:main']},
      version = "0.1.1",
      description = "Something Something bplogofun",
      long_description=long_description,
      license='GPLv3',
      url = "www.nowhere.com",
      ext_modules=[Extension('bplogofun.exact', ['src/exact.c'])],)
