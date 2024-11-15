
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cplselection

<!-- badges: start -->
<!-- badges: end -->

**cplselection** is an R package that corrects sample selection bias in
regression estimates using a Bayesian Gaussian copula method. Different
from the Heckit correction method, this method allows user to flexibly
specify the selection and the outcome marginal distributions deviating
from the bivariate normality assumption. The available options for
outcome distribution are binomial distribution (with probit and logit
link), Poisson distribution, negative binomial distribution, exponential
distribution and normal distribution.

## Installation

To install the **cplselection** package from GitHub, you need to have R
and RStudio installed on your system. Follow these steps to install the
package:

### Prerequisites

Make sure you have the `devtools` or `remotes` package installed. You
can install it from CRAN if you haven’t already:

``` r
install.packages("devtools")  # or use remotes: install.packages("remotes")
```

### Install and Load the Package

Now you are ready to install the **cplselection** package:

``` r
# Using devtools
devtools::install_github("codingcindy/cplselection")

# Alternatively, using remotes
remotes::install_github("codingcindy/cplselection")
```

After installation, you can load the package using:

``` r
library(cplselection)
```

## Usage

The **cplselection** package provides two functions for the estimation
procedure. The function `cplselectionMCMC` performs MCMC iterations and
returns an R list object containing the MCMC draws and information on
acceptance probabilities and step sizes. The other function
`cplselectionInfer` makes inferences based on the MCMC draws returned by
`cplselectionMCMC` function and returns the Bayesian estimates in an R
data.frame object.

Here’s a quick overview of how to use these functions:

### Example Usage

1.  **Generate example data set**

``` r
# Example data : Poisson outcome variable YO, logit selection process YS
library(cplselection)
nobs <- 300
dependence <- 0.5
betaO0 <- 1
betaO1 <- 1
betaS0 <- 0
betaS1 <- 1
betaS2 <- 1
uo <- runif(nobs)
us <- pnorm(dependence*qnorm(uo) + dependence*rnorm(nobs))
data <- data.frame(X1 = rnorm(nobs), X2 = rnorm(nobs))
data$YO <- qpois(p=uo, lambda=exp(betaO0 + betaO1*data$X1))
data$YS <- qbinom(p=us, prob=1/(1+exp(-(betaS0 + betaS1*data$X1 + betaS2*data$X2))), size=1)
data$YO[data$YS==0] <- NA
summary(data)
#>        X1                 X2                 YO               YS      
#>  Min.   :-3.26321   Min.   :-3.27895   Min.   : 0.000   Min.   :0.00  
#>  1st Qu.:-0.82787   1st Qu.:-0.68497   1st Qu.: 2.000   1st Qu.:0.00  
#>  Median : 0.01166   Median : 0.06364   Median : 5.000   Median :0.00  
#>  Mean   :-0.05819   Mean   : 0.01548   Mean   : 6.639   Mean   :0.49  
#>  3rd Qu.: 0.66376   3rd Qu.: 0.66409   3rd Qu.: 8.000   3rd Qu.:1.00  
#>  Max.   : 2.54044   Max.   : 2.70107   Max.   :43.000   Max.   :1.00  
#>                                        NA's   :153
```

2.  **Perform MCMC draws with `cplselectionMCMC`**

``` r
# perform MCMC draws
loop <- 10000
draws <- cplselectionMCMC(
  outcome_formula = YO~X1, select_formula = YS~X1+X2,
  outcome_dist = "Poisson", select_dist = "Logit", 
  data = data, loop = loop)
```

3.  **Bayesian inference with `cplselectionInfer`**

``` r
# if the true values are not known
result.1 <- cplselectionInfer(iterates=draws, burnin=loop/2)
print(round(result.1, 2))
#>        median   sd    lb   ub
#> rho      0.74 0.09  0.53 0.88
#> betaO0   0.94 0.07  0.79 1.09
#> betaO1   1.03 0.04  0.92 1.11
#> betaS0  -0.05 0.16 -0.38 0.26
#> betaS1   1.60 0.19  1.25 1.98
#> betaS2   1.51 0.19  1.18 1.87

# if the true values are known
trueval <- c(dependence, betaO0, betaO1, betaS0, betaS1, betaS2)
result.2 <- cplselectionInfer(iterates=draws, burnin=loop/2, trueval=trueval)
print(round(result.2, 2))
#>        trueval median   sd    lb   ub
#> rho        0.5   0.74 0.09  0.53 0.88
#> betaO0     1.0   0.94 0.07  0.79 1.09
#> betaO1     1.0   1.03 0.04  0.92 1.11
#> betaS0     0.0  -0.05 0.16 -0.38 0.26
#> betaS1     1.0   1.60 0.19  1.25 1.98
#> betaS2     1.0   1.51 0.19  1.18 1.87
```
