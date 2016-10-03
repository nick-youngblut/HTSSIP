## Test environments
* local OS X install, R 3.3.1
* ubuntu 12.04 (on travis-ci), R 3.3.1

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


### There was 1 NOTE:

~~~
.qSIP_bootstrap: no visible binding for global variable ‘OTU’
HRSIP: no visible binding for global variable ‘p’
delta_BD: no visible binding for global variable ‘SAMPLE_JOIN’
delta_BD: no visible binding for global variable ‘Count’
delta_BD: no visible binding for global variable ‘IS_CONTROL’
delta_BD: no visible binding for global variable ‘OTU’
delta_BD: no visible binding for global variable ‘Buoyant_density’
delta_BD: no visible binding for global variable ‘Count_interp’
delta_BD: no visible binding for global variable ‘center_of_mass’
... 51 lines ...
qSIP_atom_excess: no visible binding for global variable ‘IS_CONTROL’
qSIP_atom_excess: no visible binding for global variable ‘Wm’
qSIP_atom_excess: no visible binding for global variable ‘Wlab’
qSIP_atom_excess: no visible binding for global variable ‘Wlight’
qSIP_atom_excess: no visible binding for global variable ‘Gi’
qSIP_atom_excess: no visible binding for global variable ‘Mlight’
qSIP_atom_excess: no visible binding for global variable ‘Z’
qSIP_atom_excess: no visible binding for global variable ‘Mlab’
qSIP_atom_excess: no visible binding for global variable ‘Mheavymax’
qSIP_bootstrap: no visible binding for global variable ‘OTU’
qSIP_bootstrap: no visible binding for global variable ‘A’
~~~

> This note is a result of using the `%>%` declarative syntax (as described in the magrittr package) in
multiple functions in the HTSSIP package. Using this syntax improves code readability. 

