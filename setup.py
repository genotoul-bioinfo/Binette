#!/usr/bin/env python

from os import path
from setuptools import setup, find_packages
import codecs

def read(rel_path):
    here = path.abspath(path.dirname(__file__))
    with codecs.open(path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")
    

if __name__ == "__main__":
    # Get the long description from the README file
    setup_dir = path.abspath(path.dirname(__file__))
    with open(path.join(setup_dir, "README.md"), encoding="utf-8") as f:
        long_description = f.read()

    setup(
        name="binette",
        version=get_version("binette/__init__.py"),
        author="Jean Mainguy",
        packages=find_packages(),
        entry_points={"console_scripts": ["binette = binette.binette:main"]},
        url="https://github.com/genotoul-bioinfo/Binette",
        license="MIT",
        description="Binette: accurate binning refinement tool to constructs high quality MAGs.",
        long_description=(long_description),
        long_description_content_type="text/markdown",
        install_requires=[],#"pyrodigal", "pyfastx", "networkx", "checkm2"],
    )
