# ExactL3 PPI (journal version)
ExactL3 (L3E) is a better network modeling approach for link prediction in protein-protein interaction networks.

The journal publication of L3E is currently under review.

Preliminary version of L3E is published at IEEE-BIBE 2020, awarded with the Best Bioinformatics Paper Award, link:
https://doi.org/10.1109/BIBE50027.2020.00017

This repo includes the introduced algorithms, the experiments to generate all the results data, a jupyter notebook that formats all the result images, and a command-line script to run L3E and various link predictors for PPI networks.

# Requirement
Language: Python (suggested version 3.6 or above). Python Libraries (to run the program): numpy, pandas

# Usage
**Workflow: example_PPI.txt (input file) => ExactL3_cmd.py => example_PPI_ExactL3.txt (output file)**

* To run ExactL3_cmd.py in the terminal (assume Windows, cmd):
```python ExactL3_cmd.py {input file path} {output file path} {link predictor} {number of CPU core}```

* Example to run ExactL3:
```python ExactL3_cmd.py ./example_PPI.txt ./example_PPI_ExactL3.txt ExactL3 1```

* Supported Link Predictors (see the paper for reference / details): ```ExactL3```, ```L3```, ```CN```, ```CRA```, ```Sim```, ```CH2_L3```

The input file is a tab-delimited .csv file with no header, where there are two columns. The number of rows is the number of PPIs, and for each row the two columns of that row induce an non-directional PPI (each item is a protein). For example, a row ```A\tB``` implies protein 'A' and protein 'B' has a PPI.

For examples to work with our Python script, see ```./example.py```. Documentations are included as comments in the script.

# Misc
The data in our paper is generated using the script ```./src/notebook/dataGen_Yeast.py``` and ```dataGen_Human.py```, and the images in our paper are generated based on our generated data using the jupyter notebook ```./src/notebook/ppiLPred_BIBE2020_Yeast.ipynb``` and ```ppiLPred_BIBE2020_homo.ipynb```. Note that except the parsed datasets in ```./src/data/parsed/``` and generated sample PPI datasets ```./src/genData/```, no data is included since it is too large (roughly 700 GB).

# Docs
For more explanations how our python functions realize the algorithm, see [here](docs/docs.md)
