# cosme
An R Package for CoSME project


## Install package via GitHub

```r
install.packages('remotes')
remotes::install_github('bkeller2/cosme')
```

## Example script
```r
## Install package from github
# remotes::install_github('bkeller2/cosme')

## Load Package
library(cosme)

## Get Holzinger Data from `lavaan`
data(HolzingerSwineford1939, package = 'lavaan')
HSdata <- HolzingerSwineford1939[,7:15]

## Set up model via `lavaan` syntax
## Uses `lavaan::cfa()` which assumes factors
##   are correlated
m <- c(
    "visual  =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed   =~ x7 + x8 + x9"
)

## Runs 10000 reps for Information Theory
o <- cosme(
    m,           # Model Specified
    HSdata,      # Data Set
    info = TRUE  # Setting to TRUE for Info
)

## There will be 3 slots in `o`
# - `info`  information theory results
# - `freq`  frequentist results
# - `bayes` bayesian results

## Print out summaries
estimates(o)                # Estimates
estimates(o, std = TRUE)    # Standardized Estimates
fit(o)                      # Fit statistics
```

## Bug Reporting and Feature Request
Please use the [Issues](https://github.com/bkeller2/cosme/issues) tab.
