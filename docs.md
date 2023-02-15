the gold standard datasets we used are: BioGRID, STRING, MINT, HuRI. Organisms of yeast and human is supported for BioGRID, STRING, and MINT, while HuRI is a human dataset.

```python
import bioGRID as bg
bg_GI, bg_PPI = bg.parse_bioGRID(root="./src")
print(bg_PPI.head())
print(bg_GI.head())
```

for STRING, MINT, and HuRI
```python
import STRING
import MINT
import HuRI
STRING.parse_STRING(root="./src")
MINT.parse_MINT(root="./src")
HuRI.parse_HuRI(root="./src")
```

additionally we also tested the datasets in doi.org/10.1038/s41467-019-09177-y
```python
import HI_II_14_src
import IM24272_src
import Lit_BM_13_src
import Lit_NB_13_src
IM24272_src.parse_IM24272_src(root="../"),
Lit_BM_13_src.parse_Lit_BM_src(root="../"),
Lit_NB_13_src.parse_Lit_NB_src(root="../"),
HI_II_14_src.parse_HI_src(root="../")
```

The directory ```./src/notebook``` contains the jupyter notebooks that were used to generate all the data in the paper. Although not all data were uploaded onto the public space due to the size (approximately more than 5TB of effective data is generated), the sampled datasets used to generate all the precision recall curves, and the zipped data of precision recall curves, can be found in, ```./src/notebook/sampled_datasets``` and ```./src/notebook/precision_recall_out```, respectively. Below document how to unpack the data.

The basic workflow of this project is arranged as follow: sampled datasets were generated in either the notebook ```./src/notebook/generate prediction.ipynb``` or ```./src/notebook/generate prediction - human.ipynb```, and then each link predictor follow the codes in the same notebook to generate the result data into ```./src/notebook/precision_recall_out``` . The data were then processed and packed into ```./src/notebook/precision_recall_out```, then the precision recall figures generated in either ```./src/notebook/Figure - precision recall.ipynb``` or ```./src/notebook/Figure - precision recall human.ipynb```. These precision-recall notebooks included code to simply unpack the processed data to quickly regenerate the figures illustrated in the paper. The other notebooks generate data as described in their filenames.

Full methodology is described in details in - DOI: 10.1186/s12859-023-05178-3