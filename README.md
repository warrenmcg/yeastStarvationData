# yeastStarvationData

This data package contains processed data from [Marguerat et al. Cell 2012](https://doi.org/10.1016/j.cell.2012.09.019).

## Background

In this paper, fission yeast (*Schizosaccharomyces pombe*) cells were grown for 24 hours on either Edinburgh minimal media (EMM)
or EMM without NH<sub>4</sub>Cl, resulting in nitrogen starvation that leads the cells into a quiescent state. Two samples
from each group were processed for RNA-Seq using a poly-A-enriched protocol or total RNA protocol (with no poly-A-enrichment or
rRNA-depletion step). Each sample was also processed for quantification of copy numbers per cell using the NanoString assay.
Seven external controls were used to estimated copy numbers per cell of 49 mRNAs, and the results from these calibration mRNAs
were used to estimate copy numbers per cell for the rest of the yeast genome.

## Contents of the R Datasets

Contained in this data package are the following components:
+ `absCounts`: This `data.frame` contains the estimated copy numbers per cell for each of the four samples and the fold change
between the two conditions (noN vs control).

+ `denom`: This is a one-entry `character` vector with the transcript name of the yeast gene with the most consistent estimate
of copy numbers per cell across all four samples. This can be considered a validated reference gene for use in
[`sleuth-ALR`](https://github.com/warrenmcg/sleuth-ALR) or in qPCR.

+ `yeastAnnos`: This `data.frame` contains annotations for all coding and non-coding RNAs in the *S. pombe* genome, as
detailed by the [Ensembl Fungi Genomes Release 37](http://oct2017-fungi.ensembl.org/)

+ `yeastS2C`: This `data.frame` contains the metadata for the four samples in this experiment, derived from the SDRF file
from the [dataset's entry on ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1154/).

Each of the above datasets can be loaded using the following code:
```
data(absCounts, package = "yeastStarvationData")
data(denom, package = "yeastStarvationData")
data(yeastAnnos, package = "yeastStarvationData")
data(yeastS2C, package = "yeastStarvationData")
```
or alternatively
```
library(yeastStarvationData)

data(absCounts)
data(denom)
data(yeastAnnos)
data(yeastS2C)
```

## `kallisto` Datasets

This data package also contains [`kallisto`](https://github.com/pachterlab/kallisto) results that can be used for
downstream sleuth analyses. Kallisto version 0.44.0 was used with the following parameters:
```
kallisto quant -i {KALLISTO_INDEX} -b 100 -l 200 -s 20 --single \
  --fr-stranded -o {OUTPUT_FILE} {INPUT_FILE}
```

Each of the four samples has a `kallisto` directory with three components:
+ `abundance.h5`: This is an HDF5 file that `kallisto` produces, containing the expression estimates, in estimated counts and TPMs,
as well as the 100 bootstraps generated for each transcript. It is recommended to use [`sleuth`](https://github.com/pachterlab/sleuth)
or [`sleuth-ALR`](https://github.com/warrenmcg/sleuth-ALR) to load and process these files.

+ `abundance.tsv`: This is tab-separated text file containing the expression estimates (counts and TPMs) and the effective lengths
for each transcript.

+ `run_info.json`: This JSON file contains important summary statistics (e.g. % of fragments mapped) and the full command used.

These data can be accessed by the following code:
```
kal_dirs <- system.file("extdata", package = "yeastStarvationData")
## or for an individual directory
kal_dir <- system.file("extdata", "ERR135906",
  package = "yeastStarvationData", mustWork = TRUE)
```

## Using this package with [`sleuth`](https://github.com/pachterlab/sleuth)

To use this package with `sleuth`, the recommended code is the following:
```
data("yeastAnnos", package = "yeastStarvationData")
data("yeastS2C", package = "yeastStarvationData")

yeastS2C$path <- system.file("extdata", package = "yeastStarvationData")

so <- sleuth::sleuth_prep(yeastS2C, ~condition, target_mapping = yeastAnnos)
### then continue with the rest of the sleuth pipeline
```

## Using this package with [`sleuth-ALR`](https://github.com/warrenmcg/sleuth-ALR)

To use this package with `sleuth-ALR`, the recommended code is the following:
```
data("yeastAnnos", package = "yeastStarvationData")
data("yeastS2C", package = "yeastStarvationData")
data("denom", package = "yeastStarvationData")

yeastS2C$path <- system.file("extdata", package = "yeastStarvationData")

so <- sleuthALR::make_lr_sleuth_object(
        yeastS2C, ~condition,
        target_mapping = yeastAnnos,
        null_model = ~1,
        beta = "conditionnoN",
        which_var = "obs_tpm",
        denom_name = denom,
        delta = 0.01)
```

## Installation of `yeastStarvationData`

The recommended approach is to use [`devtools`](https://github.com/r-lib/devtools):
```
devtools::install_github('warrenmcg/yeastStarvationData')
```

## Full script for reproducibility

Included with this data package is a script that can rebuild all of the included datasets from scratch.
The script can be run using the following lines of code:
```
script_file <- system.file('scripts', 'make-data.R',
                 package = 'yeastStarvationData',
                 mustWork = TRUE)
source(script_file)
```

### Script Dependencies
It assumes that [`kallisto`](https://github.com/pachterlab/kallisto) is installed on the
machine's `$PATH` and that the following packages are installed:

+ `matrixStats` (`install.packages("matrixStats")`)

+ `openxlsx` (`install.packages("openxlsx")`)

### On Bioconductor
+ `Biostrings`:
```
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

The script can be run using the following line of code:
```
script_file <- system.file('scripts', 'make-data.R',
                 package = 'yeastStarvationData',
                 mustWork = TRUE)
source(script_file)
```
