# Tests 

Tests have been implemented to ensure the correctness of Binette. 


## Unit tests 

Unit tests have been implmented in the tests directory using pytest. 


To run the test suit you would need to  have install Binette from the source code.  For that, you can follow installation instructions [here](./installation.md#from-the-source-code-within-a-conda-environnement).


To install pytest in you environement you can run :

```bash
pip install .[dev]
```

Next, you can simply run the following at the root of the directory:

```bash
pytest  
```

To get the percentage of coverage of the test suit can be obtain as follow:

```bash
pytest --cov=binette 
```


```{note}

Test coverage is updated by a github workflow in the Action Tab. The test coverage report is then deployed on the github-pages and avaible [here](https://genotoul-bioinfo.github.io/Binette/). 

```


## Functional Tests


A functional test has been implemented in the CI github workflow. It performs an execution of binette on a toy dataset consisting of 4 small genomes. The test uses a checkm2 database that has been shrunk to the minimum to make diamond run faster. Finally, the results are compared with the expected results.

The test dataset is stored in this github repository: [Binette TestData](https://github.com/genotoul-bioinfo/Binette_TestData).

You can replicate this test locally by following the next steps:


1. **Install Binette**:

   Make sure you have Binette installed on your system. You can refer to the [installation](./installation.md) instructions.


2. **Clone the test dataset repository:

   Clone the dataset repository using Git:

```bash

git clone https://github.com/genotoul-bioinfo/Binette_TestData.git

cd Binette_TestData

```

3. **Run Binette**:

Run Binette on the test data with the following command

```bash
binette -b binning_results/* --contigs all_contigs.fna --checkm2_db checkm2_tiny_db/checkm2_tiny_db.dmnd -v -o test_results

```

This should complete in a few seconds.


4. **Compare Results**: 

After running Binette, you can compare the generated `final_bins_quality_reports.tsv` with the expected results stored in the `expected_results` folder. Some variation in the completeness, contamination, and score columns is expected due to Checkm2's slight variability.

You can perform the comparison manually by using the head command:

```bash
head  expected_results/final_bins_quality_reports.tsv test_results/final_bins_quality_reports.tsv

```

Alternatively, you can use the provided Python script for automated comparison: [compare_results.py](https://github.com/genotoul-bioinfo/Binette_TestData/scripts/compare_results.py) located in the scripts folder.

```bash

python scripts/compare_results.py expected_results/final_bins_quality_reports.tsv test_results/final_bins_quality_reports.tsv

```


```{warning}
The CheckM2 database used for the test dataset is very small and is only valid for the 4 genomes included in the test datasets. It should not be used elsewhere.

```