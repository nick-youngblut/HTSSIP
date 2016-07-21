HTSSIP
======
HTSSIP is an R package for analyzing high throughput 
sequence data from nucleotide stable isotope probing 
(DNA- & RNA-SIP) experiments. 

### Available analyzes 

* Identifying community-level isotope incorporatation
  * Ordinations of gradient fraction communities
  * Beta diversity of overlapping gradient fractions

* Identifying isotope incorporators
  * High resolution stable isotope probing (HR-SIP)
  * Multiple window high resolution stable isotope probing (MW-HR-SIP)
  * Quantitative stable isotope probing (q-SIP)


### Notes

* The package processes sequence data in the format of a Phyloseq object.


### Reference 

* If you use HTSSIP, please cite:


## Installation

To get the current released version from CRAN:

```R
install.packages("HTSSIP")
```

To get the current development version from github:

```R
# install.packages("devtools")
devtools::install_github("nyoungb2/HTSSIP")
```



