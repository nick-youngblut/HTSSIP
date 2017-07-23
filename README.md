HTSSIP
======

[![Travis-CI Build Status](https://travis-ci.org/nick-youngblut/HTSSIP.svg?branch=master)](https://travis-ci.org/nick-youngblut/HTSSIP)

HTSSIP is an R package for analyzing high throughput sequence data
from nucleotide stable isotope probing (DNA- & RNA-SIP) experiments. 


## Available analyzes 

* Identifying community-level isotope incorporatation
  * Ordinations of gradient fraction communities
  * Beta diversity of overlapping gradient fractions

* Identifying isotope incorporators
  * High resolution stable isotope probing (HR-SIP)
  * Multiple window high resolution stable isotope probing (MW-HR-SIP)
  * Quantitative stable isotope probing (q-SIP)


## Documentation

All documentation can be found on [CRAN](https://cran.r-project.org/package=HTSSIP).

A good place to start is the **HTSSIP introduction** vignette. 

The manuscript describing HTSSIP is:

> Youngblut, Nicholas D., Samuel E. Barnett, and Daniel H. Buckley. 2017. “HTSSIP: An R Package for Analysis of High Throughput Sequencing Data from Nucleic Acid Stable Isotope Probing (SIP) Experiments.” bioRxiv. doi:10.1101/166009.


## References 

See `References` in the **HTSSIP introduction** vignette.


## Installation

To get the current released version from [CRAN](https://cran.r-project.org/package=HTSSIP):

```R
install.packages("HTSSIP") 
```

To get the current development version from github:

```R
# install.packages("devtools")
devtools::install_github("nick-youngblut/HTSSIP")
```



