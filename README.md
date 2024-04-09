# cosme
An R Package for CoSME project


## Currently under development

This package is currently under development and has not been released yet.

## Example script
```r
## Install package from github
# remotes::install_github('bkeller2/cosme')

## Load Package
library(cosme)

## Get data
data(HolzingerSwineford1939, package = 'lavaan')
HSdata <- HolzingerSwineford1939[,7:15]

## Create Model
m <- c(
    "visual  =~ x1 + x2 + x3",
    "textual =~ x4 + x5 + x6",
    "speed   =~ x7 + x8 + x9"
)

# Runs 100 reps for Information Theory
## Set max iterations to 50
# Run
o <- cosme(
    m,
    HSdata,
    option = 
        list(
            info =  list(
                reps = 100,
                bounds = 'pos.var',
                optim.force.converged = TRUE
            )
        )
)

# There will be 3 slots in `o`
# - `info`  information theory results
# - `freq`  frequentist results
# - `bayes` bayesian results
## Then functions will be used to print out information from `o`
```

## Bug Reporting and Feature Request
Please use the [Issues](https://github.com/bkeller2/cosme/issues) tab.
