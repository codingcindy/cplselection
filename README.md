
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cplselection

<!-- badges: start -->
<!-- badges: end -->

**cplselection** is an R package that corrects sample selection bias in
regression estimates using a Bayesian Gaussian copula method. Different
from the Heckit correction method, this method allows user to flexibly
specify the selection and the outcome marginal distributions deviating
from the bivariate normality assumption.

Available options for selection distribution are

- binomial distribution (with probit link, logit link, or complementary
  log-log link)

Available options for outcome distribution are

- binomial distribution (with probit link, logit link, or complementary
  log-log link)
- Poisson distribution
- negative binomial distribution
- exponential distribution
- normal distribution

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

1.  Generate example data set

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
#>  Min.   :-2.84344   Min.   :-2.41156   Min.   : 0.000   Min.   :0.00  
#>  1st Qu.:-0.56559   1st Qu.:-0.57722   1st Qu.: 3.000   1st Qu.:0.00  
#>  Median :-0.05525   Median : 0.04793   Median : 5.000   Median :0.00  
#>  Mean   :-0.02581   Mean   : 0.03200   Mean   : 7.054   Mean   :0.49  
#>  3rd Qu.: 0.62205   3rd Qu.: 0.63440   3rd Qu.: 9.000   3rd Qu.:1.00  
#>  Max.   : 2.87384   Max.   : 2.68255   Max.   :56.000   Max.   :1.00  
#>                                        NA's   :153
```

2.  Perform MCMC draws with `cplselectionMCMC`

``` r
# Perform MCMC draws with user-specified marginal distributions
loop <- 10000
draws <- cplselectionMCMC(
  outcome_formula = YO~X1, select_formula = YS~X1+X2,
  outcome_dist = "Poisson", select_dist = "Logit", 
  data = data, loop = loop)
```

3.  Bayesian inference with `cplselectionInfer`

``` r
# If the true values are not known, returns sample selection corrected estimates:
result.1 <- cplselectionInfer(iterates=draws, burnin=loop/2)
print(round(result.1, 2))
#>        median   sd    lb   ub
#> rho      0.56 0.12  0.27 0.73
#> betaO0   1.05 0.06  0.93 1.16
#> betaO1   1.00 0.04  0.94 1.07
#> betaS0  -0.07 0.15 -0.36 0.22
#> betaS1   1.53 0.21  1.17 1.98
#> betaS2   1.26 0.19  0.93 1.68

# If the true values are known, returns true values together with sample selection 
# corrected estimates: 
trueval <- c(dependence, betaO0, betaO1, betaS0, betaS1, betaS2)
result.2 <- cplselectionInfer(iterates=draws, burnin=loop/2, trueval=trueval)
print(round(result.2, 2))
#>        trueval median   sd    lb   ub
#> rho        0.5   0.56 0.12  0.27 0.73
#> betaO0     1.0   1.05 0.06  0.93 1.16
#> betaO1     1.0   1.00 0.04  0.94 1.07
#> betaS0     0.0  -0.07 0.15 -0.36 0.22
#> betaS1     1.0   1.53 0.21  1.17 1.98
#> betaS2     1.0   1.26 0.19  0.93 1.68
```

## Documentation

More details on individual functions can be found by running:

``` r
?cplselectionMCMC

?cplselectionInfer
```

## License

This package is licensed under GPL-3
(<https://www.gnu.org/licenses/gpl-3.0.html>).
