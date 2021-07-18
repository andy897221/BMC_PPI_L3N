the gold standard datasets we used are: BioGRID, STRING, MINT, HuRI. Organisms of yeast and human is supported for BioGRID, STRING, and MINT, while HuRI is a human dataset.

```python
import bioGRID as bg
bg_GI, bg_PPI = bg.parse_bioGRID(root="./src")
print(bg_PPI.head())
print(bg_GI.head())
```

for STRING, MINT, and HuRI
```python
STRING.parse_STRING(root="./src")
MINT.parse_MINT(root="./src")
HuRI.parse_HuRI(root="./src")
```