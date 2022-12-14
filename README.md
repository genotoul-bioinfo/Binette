# Overview 



# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/GITHUB_USERNAME/bin_refinement/master/LICENSE).

# Installing

You can install bin_refinement directly from the source code or build and run it from within Docker container.

## Installing from a conda environnement

Clone this repository: 
```
git clone git@forgemia.inra.fr:jean.mainguy/bin-refinement.git
cd bin_refinement
```

Then create a Conda environment using the `bin_refinement_dev.yaml` file:
```
conda env create -n bin_refinement -f bin_refinement_dev.yaml
conda activate bin_refinement 
```

Finally install checkm2 still in bin_refinement

```
git clone --recursive https://github.com/chklovski/checkm2.git

pip install checkm2/

```
and download its database

```

checkm2 database --download --path checkm2/database/
```


Binette should be able to run :

```
python bin_refinement/bin_refinement.py -h
```


# General behaviour

## Help message

Bin_refinement can display usage information on the command line via the `-h` or `--help` argument:

```

```


# Testing

## Unit tests

You can run the unit tests for bin_refinement with the following commands:
```
cd bin_refinement/python/bin_refinement


```

## Test suite

Sample test input files are provided in the `functional_tests/test_data` folder.
```
cd functional_tests/test_data
bin_refinement two_sequence.fasta

```

Automated tests can be run using the `functional_tests/bin_refinement-test.sh` script like so:

```
cd functional_tests
./bin_refinement-test.sh -p bin_refinement -d test_data
```

The `-p` argument specifies the name of the program to test, the `-d` argument specifies the path of the directory containing test data.
The script will print the number of passed and failed test cases. More detailed information about each test case can be obtained
by requesting "verbose" output with the `-d` flag:

```
./bin_refinement-test.sh -p bin_refinement -d test_data -v
```

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[bin_refinement issue tracker](https://github.com/GITHUB_USERNAME/bin_refinement/issues)
