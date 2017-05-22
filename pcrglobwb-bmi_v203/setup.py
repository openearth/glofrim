from setuptools import setup, find_packages
import os

setup(name='pcrglobwb_bmi_v203',
      version='1.0',
      description="PCR-GLOBWB, version 2.0.3, allowing for discharge calibration at 30min resolution (only!), extended with BMI, adapted for 2-way coupling to Delft3D Flexible Mesh, converted to a library",
      long_description="""\
""",
      classifiers=[], 
      keywords='hydro',
      author='Jannis Hoch',
      author_email='j.m.hoch@uu.nl',
      url='',
      license='',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=True,
      install_requires=[
          # -*- Extra requirements: -*-
      ]
      )


