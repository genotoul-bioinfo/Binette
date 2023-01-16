#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION=''

setup(
    name='binette',
    version='0.1.0',
    author='Jean Mainguy',
    author_email='jean.mainguy@inrae.fr',
    packages=['binette'],
    package_dir={'binette': 'binette'},# entry_points={'console_scripts': ['binette = binette.binette.py']}, # url='https://github.com/GITHUB_USERNAME/binette',
    license='LICENSE',
    description=(),
    long_description=(LONG_DESCRIPTION),
    install_requires=[],
    scripts=['binette/binette.py'],
)
