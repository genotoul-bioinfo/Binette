#!/usr/bin/env python

from os import path
from setuptools import setup, find_packages

# Get the long description from the README file
setup_dir = path.abspath(path.dirname(__file__))
with open(path.join(setup_dir, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


setup(
    name="binette",
    version="0.1.1",
    author="Jean Mainguy",
    author_email="jean.mainguy@inrae.fr",
    packages=find_packages(where="src", include=["binette"]),
    package_dir={"": "src"},
    entry_points={"console_scripts": ["binette = binette.binette:main"]},
    url="https://github.com/genotoul-bioinfo/Binette",
    license="MIT",
    description="Binette: accurate binning refinement tool to constructs high quality MAGs.",
    long_description=(long_description),
    long_description_content_type="text/markdown",
    install_requires=["pyrodigal", "pyfastx"],
)
