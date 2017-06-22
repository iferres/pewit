# Pewit: Pangenome Estimation - Walks Inside Taxonomy.
This is a beta version.

## Installation
```r
library(devtools)
install_github("iferres/pewit")
```

## Dependencies

To work, it is required to have installed on your $PATH variable the following software:
 * [Mafft](http://mafft.cbrc.jp/alignment/software/)
 * [HMMER 3.1b2](http://hmmer.org/download.html)
 * [MCL](https://www.micans.org/mcl/index.html?sec_software)

### Other external dependencies

 * [Pfam-A](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/) database. The required files are: `Pfam-A.hmm.gz` and `Pfam-A.hmm.dat.gz`. Remember to decompress (gunzip) those files before runing the pipeline.

### R packages required

 * foreach
 * parallel
 * doParallel
 * utils
 * graphics
 * grDevices
 * phangorn
 * ape
 * seqinr

Some of them are part of the r-base pre-installed set of packages, so it's highly probable you have them already. The rest should be automatically downloaded and installed when installing this package.

## Note:

The PEWIT package has been designed and tested in UNIX-like platforms only.

## Citation
	"Pewit: Pangenome Estimation - Walks Inside Taxonomy"; Ignacio Ferrés, 
	Gregorio Iraola, Pablo Fresia, Daniela Costa. 
	https://github.com/iferres/pewit (2017,	Development version).
