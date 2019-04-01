#!/usr/bin/env python


from setuptools import setup, find_packages
from framework import __version__


setup(name='framework')

scripts = [
    "driver/run_meqsilhouette.py"
]


package_data = {'MeqSilhouette': [
    'input/antenna_tables/*',
    'input/source_models/*',
    'framework/*',
    'input/eht230.json',
]}


requires = [
    "numpy",
    "python_casacore",
    "meqtrees",
    "simms",
    "ATM"
    "matplotlib",
    "seaborn",
    "astropy"
]


setup(name="framework",
      version=__version__,
      description="EHT synthetic data generator",
      author="Roger Deane",
      author_email="roger.deane@up.ac.za",
      url="https://github.com/rdeane/MeqSilhouette",
      packages=find_packages(),
      package_data=package_data,
      install_requires=requires,
      scripts=scripts,
      license="limited to EHT Consortium members",
      classifiers=[],
)
