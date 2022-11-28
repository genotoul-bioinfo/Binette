# Overview 



# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/GITHUB_USERNAME/bin_refinement/master/LICENSE).

# Installing

You can install bin_refinement directly from the source code or build and run it from within Docker container.

## Installing directly from source code

Clone this repository: 
```
git clone https://github.com/GITHUB_USERNAME/bin_refinement
```

Move into the repository directory:
```
cd bin_refinement
```

Python 3 is required for this software.

Bin_refinement can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
python3 -m venv bin_refinement_dev
source bin_refinement_dev/bin/activate
pip install -U /path/to/bin_refinement
```
2. Into the global package database for all users:
```
pip install -U /path/to/bin_refinement
```
3. Into the user package database (for the current user only):
```
pip install -U --user /path/to/bin_refinement
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
