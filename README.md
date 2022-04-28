# AncestRy R package
he AncestRy package is a companion to glads (https://github.com/eriqande/glads) for spatially and temporally explicit simulation and analysis of genomic landscapes of introgression.

## Installation
To install the package run the following command in R:
```
devtools::install_github("https://github.com/LionelDiSanto/AncestRy.git", force = TRUE)
```

## Comments on the package
As currently implemented, when installing AncestRy, glads version 0.1.1 is installed (which only allow the use of the function evolve2.0, see help(evolve2.0) for details). If one would like to simulate empty demes/populations, glads version 0.1.2 needs to be installed. To that purpose, the aforementioned version of glads can be found in the *Code* directory. To install this package manually in R, run the following command line:
```
install.packages("glads_0.1.2.tar.gz", repos = NULL)
```

# Contact
For any questions or comments, contact lionel.disanto@unige.ch
