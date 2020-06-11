from setuptools import setup, find_packages
import os

install_requires = [
    'setuptools',
    'numpy', # there's an error when using 1.15.
    'bmi-python>=0.2.6', 
    'rasterio>=1.0',
    'rtree',
    'click',
    'configparser',
    'python-dateutil',
    'termcolor'
]

setup(name='glofrim',
      version='2.0',
      description="functions to couple hydrologic and hydrodynamic models using BMI",
      long_desciption="""""",
      classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development :: Libraries",
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
      ], 
      keywords='hydro',
      author='J.M. Hoch, D.M. Eilander, H. Ikeuchi',
      author_email='j.m.hoch@uu.nl',
      url='',
      license='GPLv3+',
      packages=find_packages(exclude=['docs', 'scripts', 'tests']),
      include_package_data=True,
      zip_safe=True,
      install_requires=install_requires
      )


