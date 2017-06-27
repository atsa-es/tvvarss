-   [Version](#version)
-   [Inits](#inits)
-   [MARSS model](#marss-model)
    -   [Test I: 1 observation of process](#test-i-1-observation-of-process)
    -   [Test II: 2 observations of process](#test-ii-2-observations-of-process)
    -   [Test III: 3 observations of process](#test-iii-3-observations-of-process)

### Version

This is version 0.17.06.27.

------------------------------------------------------------------------

Inits
=====

``` r
library(MARSS)
library(tvvarss)
```

    ## Loading required package: Rcpp

``` r
library(MASS)
```

MARSS model
===========

Here is a simple predator-prey model in MAR(1) form where the predator is under relatively strong density dependence, but the prey is expected to have relatively density-independent dynamics in the absence of predation.

``` r
## number of species/guilds
n_species <- 2
## number of years to simulate
n_year <- 60
## observation variance
sd_obs <- sqrt(0.05)
## initial B
BB <- matrix(0,n_species,n_species)
BB[1,1] <- 0.7
BB[1,2] <- -0.2
BB[2,1] <- 0.1
BB[2,2] <- 0.3
BB
```

    ##      [,1] [,2]
    ## [1,]  0.7 -0.2
    ## [2,]  0.1  0.3

Now we can simulate a realization of the process using `tvvarss::simTVVAR()`.

``` r
simdat <- simTVVAR(Bt = BB, topo = NULL, TT = n_year,
                   var_QX = c(0.3, 0.1), cov_QX = 0,
                   QQ_BB = matrix(0, n_species*n_species, n_species*n_species))
## other method
xx <- matrix(0, n_species, n_year)
QQ <- diag(c(0.3, 0.1))
for(t in 2:n_year) {
  xx[,t] <- BB %*% xx[,t-1] + mvrnorm(1, matrix(0, n_species, 1), QQ)
}
yy <- xx + rnorm(n_species*n_year, 0, sqrt(0.1))
```

Test I: 1 observation of process
--------------------------------

For the first test, we will use only one observation for each of the states, and set the observation variance at 0.1. We will also throw away the first 10 data points to eliminate the effect of the initial conditions.

``` r
y1 <- sim2fit(simdat, 1, sd=sd_obs, new_real = FALSE)
z1 <- t(scale(y1[,,], scale = FALSE))[,-c(1:10)]
```

Here is the model definition and fit for `MARSS()`.

``` r
## model list
mod_list <- list(
  B="unconstrained",
  U="zero",
  C="zero",
  c="zero",
  Q="diagonal and unequal",
  Z="identity",
  A="zero",
  D="zero",
  d="zero",
  R="diagonal and equal"
)
## control list
con_list <- list(maxit=9999, conv.test.slope.tol=0.1)
## fit MARSS
mf1 <- MARSS(y = z1, model = mod_list, control = con_list)
```

    ## Warning! Abstol convergence only. Maxit (=9999) reached before log-log convergence.
    ## 
    ## MARSS fit is
    ## Estimation method: kem 
    ## Convergence test: conv.test.slope.tol = 0.1, abstol = 0.001
    ## WARNING: Abstol convergence only no log-log convergence.
    ##  maxit (=9999) reached before log-log convergence.
    ##  The likelihood and params might not be at the ML values.
    ##  Try setting control$maxit higher.
    ## Log-likelihood: -55.54969 
    ## AIC: 129.0994   AICc: 131.0994   
    ##  
    ##                Estimate
    ## R.diag         0.000298
    ## B.(1,1)        0.685754
    ## B.(2,1)        0.199172
    ## B.(1,2)       -0.182662
    ## B.(2,2)        0.356133
    ## Q.(X.Y1,X.Y1)  0.248574
    ## Q.(X.Y2,X.Y2)  0.126682
    ## x0.X.Y1        0.730023
    ## x0.X.Y2       -0.787862
    ## 
    ## Standard errors have not been calculated. 
    ## Use MARSSparamCIs to compute CIs and bias estimates.
    ## 
    ## Convergence warnings
    ##  Warning: the  R.diag  parameter value has not converged.
    ##  Type MARSSinfo("convergence") for more info on this warning.

``` r
round(coef(mf1, type = "matrix")$B, 2)
```

    ##      [,1]  [,2]
    ## [1,] 0.69 -0.18
    ## [2,] 0.20  0.36

``` r
BB
```

    ##      [,1] [,2]
    ## [1,]  0.7 -0.2
    ## [2,]  0.1  0.3

Test II: 2 observations of process
----------------------------------

For the second test, we will use 2 observations for each of the states, and set the observation variance at 0.1. Again, we'll drop the first 10 observations.

``` r
y2 <- sim2fit(simdat, 2, sd=sd_obs, new_real = FALSE)
z2 <- t(scale(t(rbind(y2[,,1],y2[,,2])), scale = FALSE))[,-c(1:10)]
```

Here is the model definition and fit for `MARSS()`.

``` r
## new Z
ZZ <- matrix(0, n_species*2, n_species)
ZZ[1:2,1] <- 1
ZZ[3:4,2] <- 1
## model list
mod_list$Z <- ZZ
## fit MARSS
mf2 <- MARSS(y = z2, model = mod_list, control = con_list)
```

    ## Success! abstol and log-log tests passed at 19 iterations.
    ## 
    ## MARSS fit is
    ## Estimation method: kem 
    ## Convergence test: conv.test.slope.tol = 0.1, abstol = 0.001
    ## Estimation converged in 19 iterations. 
    ## Log-likelihood: -71.39071 
    ## AIC: 160.7814   AICc: 161.7288   
    ##  
    ##           Estimate
    ## R.diag      0.0471
    ## B.(1,1)     0.6810
    ## B.(2,1)     0.2396
    ## B.(1,2)    -0.4770
    ## B.(2,2)     0.1727
    ## Q.(X1,X1)   0.2614
    ## Q.(X2,X2)   0.0518
    ## x0.X1       0.5897
    ## x0.X2      -1.0536
    ## 
    ## Standard errors have not been calculated. 
    ## Use MARSSparamCIs to compute CIs and bias estimates.

``` r
round(coef(mf2, type = "matrix")$B, 2)
```

    ##      [,1]  [,2]
    ## [1,] 0.68 -0.48
    ## [2,] 0.24  0.17

``` r
BB
```

    ##      [,1] [,2]
    ## [1,]  0.7 -0.2
    ## [2,]  0.1  0.3

Test III: 3 observations of process
-----------------------------------

For the second test, we will use 2 observations for each of the states, and set the observation variance at 0.1. Again, we'll drop the first 10 observations.

``` r
y3 <- sim2fit(simdat, 3, sd=sd_obs, new_real = FALSE)
z3 <- t(scale(t(rbind(y3[,,1],y3[,,2])), scale = FALSE))[,-c(1:10)]
```

Here is the model definition and fit for `MARSS()`.

``` r
## new Z
ZZ <- matrix(0, n_species*3, n_species)
ZZ[1:3,1] <- 1
ZZ[4:6,2] <- 1
## model list
mod_list$Z <- ZZ
## fit MARSS
mf3 <- MARSS(y = z3, model = mod_list, control = con_list)
```

    ## Success! abstol and log-log tests passed at 16 iterations.
    ## 
    ## MARSS fit is
    ## Estimation method: kem 
    ## Convergence test: conv.test.slope.tol = 0.1, abstol = 0.001
    ## Estimation converged in 16 iterations. 
    ## Log-likelihood: -91.08432 
    ## AIC: 200.1686   AICc: 200.7893   
    ##  
    ##           Estimate
    ## R.diag      0.0488
    ## B.(1,1)     0.6431
    ## B.(2,1)     0.2208
    ## B.(1,2)    -0.2804
    ## B.(2,2)     0.1981
    ## Q.(X1,X1)   0.2948
    ## Q.(X2,X2)   0.0774
    ## x0.X1       1.2123
    ## x0.X2      -0.7387
    ## 
    ## Standard errors have not been calculated. 
    ## Use MARSSparamCIs to compute CIs and bias estimates.

``` r
round(coef(mf3, type = "matrix")$B, 2)
```

    ##      [,1]  [,2]
    ## [1,] 0.64 -0.28
    ## [2,] 0.22  0.20

``` r
BB
```

    ##      [,1] [,2]
    ## [1,]  0.7 -0.2
    ## [2,]  0.1  0.3
