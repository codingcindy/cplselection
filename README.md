
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
#>  Min.   :-3.02107   Min.   :-3.47865   Min.   : 0.000   Min.   :0.00  
#>  1st Qu.:-0.58236   1st Qu.:-0.73045   1st Qu.: 3.000   1st Qu.:0.00  
#>  Median : 0.01087   Median :-0.09963   Median : 6.000   Median :0.00  
#>  Mean   : 0.02341   Mean   :-0.04409   Mean   : 7.735   Mean   :0.49  
#>  3rd Qu.: 0.77761   3rd Qu.: 0.64504   3rd Qu.:10.000   3rd Qu.:1.00  
#>  Max.   : 3.22271   Max.   : 2.70806   Max.   :64.000   Max.   :1.00  
#>                                        NA's   :153
```

2.  Perform MCMC draws with `cplselectionMCMC`

``` r
# perform MCMC draws
loop <- 10000
draws <- cplselectionMCMC(
  outcome_formula = YO~X1, select_formula = YS~X1+X2,
  outcome_dist = "Poisson", select_dist = "Logit", 
  data = data, loop = loop)
```

3.  Bayesian inference with `cplselectionInfer`

``` r
# if the true values are not known
result.1 <- cplselectionInfer(iterates=draws, burnin=loop/2)
print(round(result.1, 2))
#>        median   sd    lb   ub
#> rho      0.69 0.10  0.45 0.83
#> betaO0   0.97 0.07  0.86 1.09
#> betaO1   1.02 0.04  0.95 1.10
#> betaS0  -0.03 0.15 -0.34 0.24
#> betaS1   1.36 0.16  1.04 1.70
#> betaS2   1.24 0.17  0.91 1.57

# if the true values are known
trueval <- c(dependence, betaO0, betaO1, betaS0, betaS1, betaS2)
result.2 <- cplselectionInfer(iterates=draws, burnin=loop/2, trueval=trueval)
print(round(result.2, 2))
#>        trueval median   sd    lb   ub
#> rho        0.5   0.69 0.10  0.45 0.83
#> betaO0     1.0   0.97 0.07  0.86 1.09
#> betaO1     1.0   1.02 0.04  0.95 1.10
#> betaS0     0.0  -0.03 0.15 -0.34 0.24
#> betaS1     1.0   1.36 0.16  1.04 1.70
#> betaS2     1.0   1.24 0.17  0.91 1.57
```

## Documentation

More details can be found by running:

``` r
?cplselectionMCMC
?cplselectionInfer
```

## License

This package is licensed under GPL-3
(<https://www.gnu.org/licenses/gpl-3.0.html>).
