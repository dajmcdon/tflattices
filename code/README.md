## Code to reproduce graphics

All code is written in `R`.

### Required packages

```r
install.packages("batchtools")
install.packages("tidyverse")
install.packages("cowplot")
install.packages("scales")
remotes::install_github("dajmcdon/tflattices", subdir = "tflatticesR")
```

### Rerunning the code

Other than the pair of `mle-v-mean-*.R` scripts, the others may be run
independently. Note that it is assumed that the home directory is set as the
root of this repository.

1. `mle-v-mean-tf-batch.R` creates the experimental results. This is intended to
   be run in parallel (or on a cluster) using the `{batchtools}` package. It
   will produce an experimental "registry" at the root, containing the results
   of the experiments.
2. `mle-v-mean-tf-results.R` collects the experimental results, does some
   processing and produces the figures in the manuscript.
