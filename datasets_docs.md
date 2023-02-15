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