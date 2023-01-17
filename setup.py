#!/usr/bin/env python

from setuptools import setup, find_packages

LONG_DESCRIPTION=''

setup(
    name='binette',
    version='0.1.1',
    author='Jean Mainguy',
    author_email='jean.mainguy@inrae.fr',
    packages=find_packages(where="src", include=['binette']), 
    package_dir={"": "src"} , # package_dir={'src': 'binette'}, 
    entry_points={'console_scripts': ['binette = binette.binette:main']}, # url='https://github.com/GITHUB_USERNAME/binette',
    license='LICENSE',
    description=(),
    long_description=(LONG_DESCRIPTION),
    install_requires=["pyrodigal", "pyfastx"], # scripts=['src/binette.py'],
)
