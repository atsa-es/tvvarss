<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/eric-ward/TVVARSS/workflows/R-CMD-check/badge.svg)](https://github.com/eric-ward/TVVARSS/actions)
[![DOI](https://zenodo.org/badge/47343421.svg)](https://zenodo.org/badge/latestdoi/47343421)
<!-- badges: end -->

<ul class="nav">

<li>

<a href="#install">Install</a>

</li>

<li>

<a href="#data">Data</a>

</li>

<li>

<a href="#models">Model structure</a>

</li>

<li>

<a href="#variances">Variances</a>

</li>

<li>

<a href="#processes">Processes</a>

</li>

<li>

<a href="#families">Families</a>

</li>

<li>

<a href="#cite">Citation info</a>

</li>

</ul>

## Installation

`tvvarss` is an R package used to simulate and/or fit time-varying
vector autoregressive state space models of community interactions.
Models are fit via Hamiltonian Monte Carlo (HMC) using the
[Stan](http://mc-stan.org/) software

After installing `devtools` (and assuming you have a C++ compiler
installed), you can install the package with

``` r
devtools::install_github("atsa-es/tvvarss")
```

A short vignette for the package is:

``` r
vignette("intro_tvvarss", package = "tvvarss")
```

### Data

`tvvarss` takes data as an array, with dimensions (site, year, species).
Using multiple sites or multiple species is completely optional, and we
assume 1 measurement per site-year-species combination (for multiple
measurements, these may be created as replicate sites).

Simulated data may be generated with our built in simulation function
`simTVVAR`, which is demonstrate in the vignette

### Model structure

We allow the \(\textbf{B}\) matrix in a conventional MAR model to be
either static or time varying – time varying models are assumed by
default, but may changed to static matrices with the

``` r
tvvarss(y, dynamicB = FALSE)
```

Because the number of elements in \(\textbf{B}\) may not be
identifiable, we also allow constraints to be passed in, depending on
the food web structure. This is done with the `topo` argument, which is
detailed in the package vignette.

### Variances

We allow process and observation variances to be shared across sites or
species. These are set with the arguments is dimensioned (number of
species x number of sites), and elements are passed in as integers.

``` r
tvvarss(y, shared_q = Q, shared_r = R)
```

For example, if we had a dataset consiting of 3 sites and 5 species, the
default form of these matrices in the function would be a unique
variance by species, that was shared across sites, e.g.

    m = matrix(0, 5, 3)
    for(i in 1:5) m[i,] = i
    colnames(m) = paste0("Site ",1:3)
    rownames(m) = paste0("Species ", 1:5)
    knitr::kable(m)

We could change this up, and make variances be unique to each
species-site

    m = matrix(1:15, 5, 3)
    colnames(m) = paste0("Site ",1:3)
    rownames(m) = paste0("Species ", 1:5)
    knitr::kable(m)

Or equal across species - sites

    m = matrix(1, 5, 3)
    colnames(m) = paste0("Site ",1:3)
    rownames(m) = paste0("Species ", 1:5)
    knitr::kable(m)

### Processes

In some situations, like the replicate observations described above, we
may want to relate sites as measurements of the same underlying state.
This can be done with the ‘process’ argument, which maps sites to
states.

``` r
tvvarss(y = y, process = Z)
```

By default each site is its own process, so for the 3 site example,
`process` would be

``` r
process = 1:3
```

These could be shared though – if sites 1 and 3 were a different process
than site 2, this could be expressed as

``` r
process = c(1,2,1)
```

### Families

We’ve implemented a few distribution families that may be useful for
modeling non-normal data.

The default family is “gaussian”, but

``` r
tvvarss(y = y, family="gaussian")
```

the family can be any of the following: “gaussian”, “binomial”,
“poisson”, “gamma”, “lognormal”. We assumea logit link for modeling
the binomial, and a log link for the Gamma, poisson, and lognormal.

### CITATION

Ward, E.J., M.D. Scheuerell, and S.L. Katz. 2021. ‘tvvarss’: Time
varying vector autoregressive models in Stan.
[![DOI](https://zenodo.org/badge/47343421.svg)](https://zenodo.org/badge/latestdoi/47343421)
