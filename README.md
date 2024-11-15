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

    install.packages("devtools")  # or use remotes: install.packages("remotes")

### Install and Load the Package

Now you are ready to install the **cplselection** package:

    # Using devtools
    devtools::install_github("codingcindy/cplselection")

    # Alternatively, using remotes
    remotes::install_github("codingcindy/cplselection")

After installation, you can load the package using:

    library(cplselection)

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

1.  Perform MCMC draws with `cplselectionMCMC`

<!-- -->

    # Example data : Poisson outcome variable YO, logit selection process YS
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

    # perform MCMC draws
    loop <- 10000
    draws <- cplselectionMCMC(
      outcome_formula = YO~X1, select_formula = YS~X1+X2,
      outcome_dist = "Poisson", select_dist = "Logit", 
      data = data, loop = loop)
    #>     sampling progress [>--------------------------------]   2% in  0s | eta: 10s    sampling progress [>--------------------------------]   3% in  0s | eta: 10s    sampling progress [>--------------------------------]   3% in  0s | eta: 11s    sampling progress [>--------------------------------]   4% in  0s | eta: 11s    sampling progress [>--------------------------------]   4% in  0s | eta: 10s    sampling progress [>--------------------------------]   4% in  0s | eta: 11s    sampling progress [>--------------------------------]   4% in  0s | eta: 10s    sampling progress [>--------------------------------]   4% in  0s | eta: 11s    sampling progress [>--------------------------------]   4% in  0s | eta: 10s    sampling progress [>--------------------------------]   4% in  0s | eta: 11s    sampling progress [>--------------------------------]   4% in  0s | eta: 10s    sampling progress [>--------------------------------]   4% in  0s | eta: 11s    sampling progress [>--------------------------------]   4% in  1s | eta: 11s    sampling progress [>--------------------------------]   5% in  1s | eta: 11s    sampling progress [=>-------------------------------]   5% in  1s | eta: 11s    sampling progress [=>-------------------------------]   6% in  1s | eta: 11s    sampling progress [=>-------------------------------]   7% in  1s | eta: 11s    sampling progress [=>-------------------------------]   8% in  1s | eta: 11s    sampling progress [==>------------------------------]   8% in  1s | eta: 11s    sampling progress [==>------------------------------]   9% in  1s | eta: 11s    sampling progress [==>------------------------------]  10% in  1s | eta: 11s    sampling progress [==>------------------------------]  10% in  1s | eta: 10s    sampling progress [==>------------------------------]  10% in  1s | eta: 11s    sampling progress [==>------------------------------]  10% in  1s | eta: 10s    sampling progress [==>------------------------------]  11% in  1s | eta: 10s    sampling progress [===>-----------------------------]  11% in  1s | eta: 10s    sampling progress [===>-----------------------------]  12% in  1s | eta: 10s    sampling progress [===>-----------------------------]  13% in  1s | eta: 10s    sampling progress [===>-----------------------------]  13% in  2s | eta: 10s    sampling progress [===>-----------------------------]  14% in  2s | eta: 10s    sampling progress [====>----------------------------]  14% in  2s | eta: 10s    sampling progress [====>----------------------------]  15% in  2s | eta: 10s    sampling progress [====>----------------------------]  16% in  2s | eta: 10s    sampling progress [====>----------------------------]  17% in  2s | eta: 10s    sampling progress [=====>---------------------------]  17% in  2s | eta: 10s    sampling progress [=====>---------------------------]  18% in  2s | eta: 10s    sampling progress [=====>---------------------------]  19% in  2s | eta: 10s    sampling progress [=====>---------------------------]  19% in  2s | eta:  9s    sampling progress [=====>---------------------------]  20% in  2s | eta:  9s    sampling progress [======>--------------------------]  20% in  2s | eta:  9s    sampling progress [======>--------------------------]  21% in  2s | eta:  9s    sampling progress [======>--------------------------]  21% in  3s | eta:  9s    sampling progress [======>--------------------------]  22% in  3s | eta:  9s    sampling progress [======>--------------------------]  23% in  3s | eta:  9s    sampling progress [=======>-------------------------]  23% in  3s | eta:  9s    sampling progress [=======>-------------------------]  24% in  3s | eta:  9s    sampling progress [=======>-------------------------]  25% in  3s | eta:  9s    sampling progress [=======>-------------------------]  26% in  3s | eta:  9s    sampling progress [========>------------------------]  26% in  3s | eta:  9s    sampling progress [========>------------------------]  27% in  3s | eta:  9s    sampling progress [========>------------------------]  28% in  3s | eta:  9s    sampling progress [========>------------------------]  28% in  3s | eta:  8s    sampling progress [========>------------------------]  29% in  3s | eta:  8s    sampling progress [=========>-----------------------]  29% in  3s | eta:  8s    sampling progress [=========>-----------------------]  30% in  3s | eta:  8s    sampling progress [=========>-----------------------]  30% in  4s | eta:  8s    sampling progress [=========>-----------------------]  31% in  4s | eta:  8s    sampling progress [=========>-----------------------]  32% in  4s | eta:  8s    sampling progress [==========>----------------------]  32% in  4s | eta:  8s    sampling progress [==========>----------------------]  33% in  4s | eta:  8s    sampling progress [==========>----------------------]  34% in  4s | eta:  8s    sampling progress [==========>----------------------]  35% in  4s | eta:  8s    sampling progress [===========>---------------------]  35% in  4s | eta:  8s    sampling progress [===========>---------------------]  36% in  4s | eta:  8s    sampling progress [===========>---------------------]  36% in  4s | eta:  7s    sampling progress [===========>---------------------]  37% in  4s | eta:  7s    sampling progress [===========>---------------------]  38% in  4s | eta:  7s    sampling progress [============>--------------------]  38% in  4s | eta:  7s    sampling progress [============>--------------------]  38% in  5s | eta:  7s    sampling progress [============>--------------------]  38% in  5s | eta:  8s    sampling progress [============>--------------------]  38% in  5s | eta:  7s    sampling progress [============>--------------------]  39% in  5s | eta:  7s    sampling progress [============>--------------------]  40% in  5s | eta:  7s    sampling progress [============>--------------------]  41% in  5s | eta:  7s    sampling progress [=============>-------------------]  41% in  5s | eta:  7s    sampling progress [=============>-------------------]  42% in  5s | eta:  7s    sampling progress [=============>-------------------]  43% in  5s | eta:  7s    sampling progress [=============>-------------------]  44% in  5s | eta:  7s    sampling progress [==============>------------------]  44% in  5s | eta:  7s    sampling progress [==============>------------------]  45% in  5s | eta:  7s    sampling progress [==============>------------------]  46% in  5s | eta:  7s    sampling progress [==============>------------------]  46% in  6s | eta:  7s    sampling progress [==============>------------------]  46% in  6s | eta:  6s    sampling progress [==============>------------------]  47% in  6s | eta:  6s    sampling progress [===============>-----------------]  47% in  6s | eta:  6s    sampling progress [===============>-----------------]  48% in  6s | eta:  6s    sampling progress [===============>-----------------]  49% in  6s | eta:  6s    sampling progress [===============>-----------------]  50% in  6s | eta:  6s    sampling progress [================>----------------]  50% in  6s | eta:  6s    sampling progress [================>----------------]  51% in  6s | eta:  6s    sampling progress [================>----------------]  52% in  6s | eta:  6s    sampling progress [================>----------------]  53% in  6s | eta:  6s    sampling progress [=================>---------------]  53% in  6s | eta:  6s    sampling progress [=================>---------------]  54% in  6s | eta:  6s    sampling progress [=================>---------------]  54% in  6s | eta:  5s    sampling progress [=================>---------------]  54% in  7s | eta:  5s    sampling progress [=================>---------------]  55% in  7s | eta:  5s    sampling progress [=================>---------------]  56% in  7s | eta:  5s    sampling progress [==================>--------------]  56% in  7s | eta:  5s    sampling progress [==================>--------------]  57% in  7s | eta:  5s    sampling progress [==================>--------------]  58% in  7s | eta:  5s    sampling progress [==================>--------------]  59% in  7s | eta:  5s    sampling progress [===================>-------------]  59% in  7s | eta:  5s    sampling progress [===================>-------------]  60% in  7s | eta:  5s    sampling progress [===================>-------------]  61% in  7s | eta:  5s    sampling progress [===================>-------------]  62% in  7s | eta:  5s    sampling progress [====================>------------]  62% in  7s | eta:  5s    sampling progress [====================>------------]  62% in  7s | eta:  4s    sampling progress [====================>------------]  63% in  7s | eta:  4s    sampling progress [====================>------------]  63% in  8s | eta:  4s    sampling progress [====================>------------]  64% in  8s | eta:  4s    sampling progress [====================>------------]  65% in  8s | eta:  4s    sampling progress [=====================>-----------]  65% in  8s | eta:  4s    sampling progress [=====================>-----------]  66% in  8s | eta:  4s    sampling progress [=====================>-----------]  67% in  8s | eta:  4s    sampling progress [=====================>-----------]  68% in  8s | eta:  4s    sampling progress [======================>----------]  68% in  8s | eta:  4s    sampling progress [======================>----------]  69% in  8s | eta:  4s    sampling progress [======================>----------]  70% in  8s | eta:  4s    sampling progress [======================>----------]  70% in  8s | eta:  3s    sampling progress [======================>----------]  71% in  8s | eta:  3s    sampling progress [=======================>---------]  71% in  8s | eta:  3s    sampling progress [=======================>---------]  72% in  8s | eta:  3s    sampling progress [=======================>---------]  72% in  9s | eta:  3s    sampling progress [=======================>---------]  73% in  9s | eta:  3s    sampling progress [=======================>---------]  74% in  9s | eta:  3s    sampling progress [========================>--------]  74% in  9s | eta:  3s    sampling progress [========================>--------]  75% in  9s | eta:  3s    sampling progress [========================>--------]  76% in  9s | eta:  3s    sampling progress [========================>--------]  77% in  9s | eta:  3s    sampling progress [=========================>-------]  77% in  9s | eta:  3s    sampling progress [=========================>-------]  78% in  9s | eta:  3s    sampling progress [=========================>-------]  79% in  9s | eta:  3s    sampling progress [=========================>-------]  79% in  9s | eta:  2s    sampling progress [=========================>-------]  80% in  9s | eta:  2s    sampling progress [==========================>------]  80% in  9s | eta:  2s    sampling progress [==========================>------]  80% in 10s | eta:  2s    sampling progress [==========================>------]  81% in 10s | eta:  2s    sampling progress [==========================>------]  82% in 10s | eta:  2s    sampling progress [==========================>------]  83% in 10s | eta:  2s    sampling progress [===========================>-----]  83% in 10s | eta:  2s    sampling progress [===========================>-----]  84% in 10s | eta:  2s    sampling progress [===========================>-----]  85% in 10s | eta:  2s    sampling progress [===========================>-----]  86% in 10s | eta:  2s    sampling progress [============================>----]  86% in 10s | eta:  2s    sampling progress [============================>----]  87% in 10s | eta:  2s    sampling progress [============================>----]  87% in 10s | eta:  1s    sampling progress [============================>----]  88% in 10s | eta:  1s    sampling progress [============================>----]  89% in 10s | eta:  1s    sampling progress [============================>----]  89% in 11s | eta:  1s    sampling progress [=============================>---]  89% in 11s | eta:  1s    sampling progress [=============================>---]  90% in 11s | eta:  1s    sampling progress [=============================>---]  91% in 11s | eta:  1s    sampling progress [=============================>---]  92% in 11s | eta:  1s    sampling progress [==============================>--]  92% in 11s | eta:  1s    sampling progress [==============================>--]  93% in 11s | eta:  1s    sampling progress [==============================>--]  94% in 11s | eta:  1s    sampling progress [==============================>--]  95% in 11s | eta:  1s    sampling progress [===============================>-]  95% in 11s | eta:  1s    sampling progress [===============================>-]  96% in 11s | eta:  1s    sampling progress [===============================>-]  96% in 11s | eta:  0s    sampling progress [===============================>-]  97% in 11s | eta:  0s    sampling progress [===============================>-]  98% in 11s | eta:  0s    sampling progress [===============================>-]  98% in 12s | eta:  0s    sampling progress [================================>]  98% in 12s | eta:  0s    sampling progress [================================>]  99% in 12s | eta:  0s    sampling progress [================================>] 100% in 12s | eta:  0s    sampling progress [=================================] 100% in 12s | eta:  0s

1.  Bayesian inference with `cplselectionInfer`

<!-- -->

    # if the true values are not known
    result.1 <- cplselectionInfer(iterates=draws, burnin=loop/2)
    print(round(result.1, 2))
    #>        median   sd    lb   ub
    #> rho      0.73 0.08  0.55 0.84
    #> betaO0   0.94 0.05  0.86 1.03
    #> betaO1   1.02 0.03  0.95 1.06
    #> betaS0  -0.14 0.14 -0.40 0.15
    #> betaS1   1.51 0.19  1.14 1.90
    #> betaS2   1.24 0.16  0.94 1.57

    # if the true values are known
    trueval <- c(dependence, betaO0, betaO1, betaS0, betaS1, betaS2)
    result.2 <- cplselectionInfer(iterates=draws, burnin=loop/2, trueval=trueval)
    print(round(result.2, 2))
    #>        trueval median   sd    lb   ub
    #> rho        0.5   0.73 0.08  0.55 0.84
    #> betaO0     1.0   0.94 0.05  0.86 1.03
    #> betaO1     1.0   1.02 0.03  0.95 1.06
    #> betaS0     0.0  -0.14 0.14 -0.40 0.15
    #> betaS1     1.0   1.51 0.19  1.14 1.90
    #> betaS2     1.0   1.24 0.16  0.94 1.57
