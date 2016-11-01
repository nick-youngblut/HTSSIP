## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1
* win-builder (devel)

## R CMD check results

### There were no ERRORS

### There was 1 WARNING:

~~~
Warning: replacing previous import ‘BiocGenerics::Position’ by ‘ggplot2::Position’ when loading ‘phyloseq’
All user-level objects in a package should have documentation entries.
See chapter ‘Writing R documentation files’ in the ‘Writing R
Extensions’ manual.
~~~

> The warning is a result of including the phyloseq R package in the list of Imports (see DESCRIPTION).


### There was 1 NOTE from win-builder:

Note 1:

Maintainer: ‘Nicholas Youngblut <nyoungb2@gmail.com>’
New submission


## Other comments

* The Windows build error: "cannot open file 'startup.Rs': No such file or directory" should be fixed.
* We have reduced the run times for both building vignettes and running tests. 
* We have made the package description more descriptive and included references with DOIs.

