# ExactL3 Link Prediction for PPI Networks
ExactL3 (L3E) is a better network modeling approach for link prediction in protein-protein interaction networks.

The journal publication of L3E is currently under review.

Preliminary version of L3E is published at IEEE-BIBE 2020, awarded with the Best Bioinformatics Paper Award, link:
https://doi.org/10.1109/BIBE50027.2020.00017

This repo includes the proposed algorithms, jupyter notebooks that code all the experiments for data & figure generations, and a command-line script to run L3 and L3E link predictors for PPI networks.

# Requirement
Language: ```python``` (suggested version 3.6 or above).

Python libraries requirement: ```numpy```, ```pandas```

# Usage
To run main.py in the terminal:

```python main.py {input file path} {output file path} {link predictor}```

Example to run L3E(f1):

```python main.py ./sample_data/example_PPI.txt ./sample_data/out.txt L3EPrime(f1)```

Supported Link Predictors (see the paper details): ```L3```, ```L3EPrime(f1)```, ```L3EPrime(f2)```, ```L3E(f1)```, ```L3E(f2)```

The input file is a tab-delimited .csv file with no header & with two columns. The number of rows is the amount of PPIs. Each row represents a pair of protein (two columns) that exists a PPI. For example, a row ```A\tB``` implies protein 'A' and protein 'B' has a PPI.

To use L3E programmatically, see ```./example.py``` as an example. It is encouraged to use L3E this way for real PPI datasets, so thata custom multiprocessing handler can be coded and parsed into the L3E function.

# Implementations
```core.py``` is the actual implementation of the L3E and L3 algorithm. The code is self-explanatory.

```example.py``` and ```main.py```: see [Usage](#Usage).

```unitTest.py``` provides basic unit test of simple types of PPI graphs. it is to show that the predictors implementation in ```core.py``` is the same as the predictors in ```./src/PPILinkPred.py```, of which ```./src/PPILinkPred.py``` is the actual implementation used for the experiments in the paper (e.g. multiprocessing handler is included to generate data).

# Documentations
Details of L3E networking modeling is available in the paper (under review), and details of the folder ```./src/notebook``` (data processing, results, and figures) are elaborated [here](docs/docs.md)