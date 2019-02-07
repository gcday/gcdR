# gcdR
Useful wrappers and data structures for simple and reproducible processing of Seurat objects, from count-matrices all the way to downstream differential expression and gene set enrichment analyses. 


## Installation
Important: restart R before running these commands! (under the Session tab in RStudio)
I believe this requires R 3.5 or greater! If this process doesn't work, I reccomend downloading the latest version [here](http://cran.cnr.berkeley.edu/). 
```{r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("satijalab/seurat", ref = "release/3.0", force = T)
devtools::install_github("gcday/gcdR")
```
Note -- to calculate UMAP reduction in addition to t-SNE, you'll need to install UMAP by running the following command in Terminal (on MacOS) or Command Prompt in Windows:
```
pip install umap-learn
```
The above command relies on Python being installed. If this command fails, one reason may be that Python isn't installed. (You can check by running ```python``` on the command line). Try installing the latest Python [here](https://www.python.org/downloads/release/python-372/)

The Shiny app will work just fine without the ```umap-learn``` package, but you'll be limited to t-SNE. 

## Quick run
```{r}
library(gcdR)
runCellViewer()
```



