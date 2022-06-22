##packaging command python setup.py sdist

import sys
import os.path
from setuptools import setup, Extension, find_packages


setup(
    name='Locateguide',
    version='0.3',
    packages=find_packages(),
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.txt').read(),
    package_dir = {'Locateguide': 'Locateguide'},
    package_data={'test': ['Locateguide/test/*']},
    py_modules = ['Locateguide.Detect_peaks', 'Locateguide.full_analysis', 'Locateguide.Validate'],
    install_requires=['numpy==1.22.0' ,
                      'pyfaidx==0.5.1', 
                      'subprocess32==3.2.7' ,
                        'swalign==0.3.4',
                        'HTSeq==0.9.1',
                        'statsmodels==0.8.0']

)


