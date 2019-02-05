# gcdR
Useful wrappers and data structures for simple and reproducible processing of Seurat objects, from count-matrices all the way to downstream differential expression and gene set enrichment analyses. 


## Installation
```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("satijalab/seurat", ref = "release/3.0", force = T)
devtools::install_github("gcday/gcdR")
```
Note -- to calculate UMAP reduction in addition to t-SNE, you'll need to install UMAP by running the following command in Terminal (on MacOS) or Command Prompt in Windows:
```
pip install umap-learn
```
Note -- this requires Python to be installed. If this command fails, try installing the latest Python [here](https://www.python.org/downloads/release/python-372/)


## Quick run
```{r}
library(gcdR)
runCellViewer()
```



