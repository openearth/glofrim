from setuptools import setup, find_packages
import os

setup(name='glofrim',
      version='1.1',
      description="functions to couple hydrology with 1D hydrodynamics",
      long_desciption="""\
""",
      classifiers=[], 
      keywords='hydro',
      author='J.M. Hoch, D.M. Eilander, H. Ikeuchi',
      author_email='jannis.hoch@deltares.nl',
      url='',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ]
      )


