
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MMRFVariantPackage

<!-- badges: start -->

<!-- badges: end -->

The goal of MMRFVariantPackage is to analize variants from MMRF-CoMMpass
dataset.

## Installation

Once R (version \> “4.0”) has been started, you can install the released
version of MMRFVariantPackage from GitHub with:

``` r
devtools::install_github("marziasettino/MMRFVariantPackage", build_vignettes = TRUE)
library(MMRFVariantPackage)
```

## Required libraries

``` r
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)
```

## Vignettes

A list of all currently integrated vignettes can be obtained through:

``` r
vignette(package="MMRFVariant")
```

The best way to view vignettes is in your web browser:

``` r
devtools::load_all(".")
browseVignettes("MMRFVariant")
```

Get the list of the example data sets

``` r
data(package = "MMRFVariant")
```
