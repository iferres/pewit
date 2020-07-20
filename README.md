  <!-- badges: start -->
  [![R build status](https://github.com/iferres/pewit/workflows/R-CMD-check/badge.svg)](https://github.com/iferres/pewit/actions)
  [![codecov](https://codecov.io/gh/iferres/pewit/branch/master/graph/badge.svg)](https://codecov.io/gh/iferres/pewit)
  <!-- badges: end -->

# Pewit: Pangenome Estimation - Walks Inside Taxonomy.
This is a beta version.

## Installation
```r
library(devtools)
install_github("iferres/pewit")
```

## Dependencies

To work, it is required to have installed on your $PATH variable the following software:
 * [HMMER 3.1b2](http://hmmer.org/download.html)
 * [MCL](https://www.micans.org/mcl/index.html?sec_software)

### Other external dependencies

 * [Pfam-A](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/) database. The required files are: `Pfam-A.hmm.gz` and `Pfam-A.hmm.dat.gz`. Remember to decompress (gunzip) those files before running the pipeline.


## Note:

The PEWIT package has been designed and tested in UNIX-like platforms only.

## Citation
	"Pewit: Pangenome Estimation - Walks Inside Taxonomy"; Ignacio Ferres, Gregorio Iraola, Pablo Fresia, Daniela Costa. 
	https://github.com/iferres/pewit (2017,	Development version).
