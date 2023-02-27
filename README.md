# Normalized L3-based Link Prediction in Protein-Protein Interaction Networks
NormalizedL3 (L3N) is a network modeling approach for link prediction in protein-protein interaction networks.

The journal publication of L3N has currently online in BMC Bioinformatics: [10.1186/s12859-023-05178-3](https://doi.org/10.1186/s12859-023-05178-3)

This repo includes the proposed algorithms, jupyter notebooks that code all the experiments for data & figure generations, and a command-line script to run L3 and L3N link predictors for PPI networks.

# Requirement
Language: ```python``` (version 3.6 or above).

Python libraries requirement: ```numpy```, ```pandas```

# Usage
To run main.py in the terminal:

```python main.py {input file path} {output file path} {link predictor}```

Example to run L3N(f1):

```python main.py ./sample_data/example_PPI.txt ./sample_data/out.txt L3NPrime(f1)```

Supported Link Predictors (see the paper details): ```L3```, ```L3NPrime(f1)```, ```L3NPrime(f2)```, ```L3N(f1)```, ```L3N(f2)```

The input file is a tab-delimited .csv file with no header & with two columns. The number of rows is the amount of PPIs. Each row represents a pair of protein (two columns) that exists a PPI. For example, a row ```A\tB``` implies protein 'A' and protein 'B' has a PPI.

To use L3N programmatically, see ```./example.py``` as an example. It is encouraged to use L3E this way for real PPI datasets, so that custom multiprocessing handler can be coded and parsed into the L3E function.

# Implementations
```core.py``` is the actual implementation of the L3N and L3 algorithm. The code is self-explanatory.

```example.py``` and ```main.py```: see [Usage](#Usage).

```unitTest.py``` provides basic unit tests of simple types of PPI graphs. it is to show that the predictors implementation in ```core.py``` is the same as the predictors in ```./src/PPILinkPred.py``` under the unit tests, of which ```./src/PPILinkPred.py``` is the actual implementation used for the experiments in the paper (e.g. multiprocessing handler is included to generate data).

# Documentations
Details to use the datasets that generated the data in our publication is documented [here](datasets_docs.md)

Details of L3N networking modeling is available in the paper, [10.1186/s12859-023-05178-3](https://doi.org/10.1186/s12859-023-05178-3)

Details of the folder ```./src/notebook``` (data processing, results, and data figures) are briefly elaborated [here](docs.md)
