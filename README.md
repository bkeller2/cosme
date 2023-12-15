# cosme
An R Package for CoSME project


## Currently under development

This package is currently under development and has not been released yet.

## Example script
```r

## Install and load package
remotes::install_github('bkeller2/cosme')
library(cosme)

## Place holder to deal with model
cosme_model <- function(model) {
    # TODO check if constructed appropriately
    
    # Return object
    structure(model, class = 'cosme_model')
}

## Get data
data(HolzingerSwineford1939, package = 'lavaan')
HSdata <- HolzingerSwineford1939[,7:15]

## Create Model
m <- cosme_model(
    "visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9"
)

# Set up multisession evaluation
# Should speed up computations in general
future::plan('multisession')

# Runs 25 reps for testing with ockhamSEM
## Additional functional inputs will be added for `cosme`
## To set up models for each method
o <- cosme(m, HSdata)

# There will be 3 slots in `o`
# - `info`  information theory results
# - `freq`  frequentist results
# - `bayes` bayesian results
## Then functions will be used to print out information from `o`
```

## Bug Reporting and Feature Request
Please use the [Issues](https://github.com/bkeller2/cosme/issues) tab.
