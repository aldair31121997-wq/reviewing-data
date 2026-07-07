# HPAIVpredict

`HPAIVpredict` is designed to predict the aquisition of insertions in a given influenza H5/H7 RNA sequence, according to our slippage model. Two steps are required:

* the first is the calculation of backtrack DG values of each position (using a python script);
* the second is the use of a multivariate regresion model on the first step results, to obtain insertions that are likely to occur and an associated graphical representation.


## Dependencies

To run `HPAIVpredict`, please ensure that the following dependencies are satisfied.

## Python packages

The script `inserpredictor.py` script was tested on Python 3.7 but should be fine with more recent versions of Python.

* [ViennaRNA](https://pypi.org/project/ViennaRNA/) (version 2.7.2)
* [pandas](https://pandas.pydata.org/) (version 3.0.3)
* [numpy](https://numpy.org/) (version 2.4.4)
* [statsmodels](https://www.statsmodels.org/stable/api.html) (version 0.14.6)
* [Biopython module Bio.Seq](https://biopython.org/docs/1.75/api/Bio.Seq.html) (version 1.74)

an updated requirements file is aviable at requirements/requirements.txt, to use it run 

```
pip install -r requirements/requirements.txt

```


### R packages

* [readxl](https://readxl.tidyverse.org/) (version 1.5.0)
* [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html) (version 4.6.0)
* [ggplot2](https://ggplot2.tidyverse.org/) (version 4.0.3)
* [ggrepel](https://ggrepel.slowkow.com/) (version 0.9.8)
* [ggthemes](https://ggplot2.tidyverse.org/reference/ggtheme.html) (version 4.0.3)
* [dplyr](https://dplyr.tidyverse.org/) (version 1.2.1)
* [bioseq](https://cran.r-project.org/web/packages/bioseq/index.html) (version 0.1.5)
* [stringi](https://cran.r-project.org/web/packages/stringi/index.html) (version 1.8.7)
* [writexl](https://cran.r-project.org/web/packages/writexl/index.html) (version 1.5.4)
* [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html) (version 0.6.3)
* [rstudio](https://docs.posit.co/ide/user/) (version 2023.12.1.402)


## 1st step: Get backtrack scores

The sequence used for the prediction must be available in `fasta`. It is recommended that only small regions of a given RNA sequence are predicted since the algorithm grows exponentially with the sequence length and the number of sequences to predict. The step is performed using:

```
inserpredictor.py PATH_TO_SEQUENCE_FILE.fasta start stop correction
```

with:

| **argument** | **value** |
| ------------ | --------- |
| `start`      | sequence position at which to start the algorithm |
| `stop`       | sequence position at which to stop the algorithm |
| `correction` | new position assigned to the start position |


**Example**:

```
python3 inserpredictor.py data/allviruses.fasta 26 85 986
```

**WHERE TO FIND allviruses.fasta**

the HA partial sequences used on this study, can be found at the data/ folder under the name allviruses.fasta

**Value**:

The script creates two files (`crnaDDG.xlsx` and `vrnaDDG.xlsx`). The first corresponds to the backtrack positions during vRNA syntesis (with cRNA as the template), while the second corresponds to the backtrack positions during cRNA syntesis (with vRNA as the template).


## 2nd step: Get prediction frequencies and scores

We recommend to use [RStudio](https://posit.co/products/open-source/rstudio) to run the file `visualizator.Rmd`, the only thing you need to modify are the following lines:

```
predictionH5 <- "PATH_TO_crnaDDG.xlsx"
predictionH5crna <- "PATH_TO_vrnaDDG.xlsx"
```

where these two files are those obtained in the first step. The script outputs a table in `xlsx` format containing all the analyzed positions and estimated likelihoods of insertions. When a given insertion has a prediction score higher than the detection threshold ($10^{-5}$), the prediction column contains its estimated $\log_{10}$-frequency unless the position is filtered by the `slippage accesibilty` rule. In cases where the prediction score is lower than the detection threshold or is filtered by the `slippage accesibilty` rule, the prediction column contains `BELOWtreshold`. 

In addition, the script creates the same type of plot than the plots shown in the article using a given window of visualization. this window corresponds to the positions in the HA 
predicted and obseved (going from nucleotide position 1012 to 1030)


## Additional information on the model fitting

Model training and fitting are available in the `training-fitting/` folder, which contains the R script used to generate the linear regression model, as well as the table containing the insertions and their corresponding ΔΔG, ΔG, and match/mismatch values used for training.


# Note

**TODO: add preprint citation**

This version of `HPAIVpredict` (v 1.0.0) is the first to be published but not the last: A new version including substitution forecasting and furin-clevage scores is currently being developed. It will be released as a single user-friendly module and available via this repository send a mail to [aldair31121997@gmail.com](mailto:aldair3112
