# leiv: Bivariate Linear Errors-In-Variables Estimation
***
This repository contains SAS macros for estimating the slope and intercept of a bivariate linear relationship by calculating a posterior density that is invariant to interchange and scaling of the coordinates.

The SAS macros implement the method described in [Leonard, David. Estimating a bivariate linear relationship. Bayesian Anal. 6 (2011), no. 4, 727--754. doi:10.1214/11-BA627. https://projecteuclid.org/euclid.ba/1339616542](https://projecteuclid.org/euclid.ba/1339616542) and are a companion to the R package, [https://CRAN.R-project.org/package=leiv](https://CRAN.R-project.org/package=leiv).

## Contents
- [leiv_macros.sas](leiv_macros.sas) (SAS macros for calculating the _leiv_ posterior density and related estimates)
- [leiv_examples.sas](leiv_examples.sas) (SAS program illustrating calls to the _leiv_ SAS macro.)

## Requirements
- Base SAS
- SAS/STAT, version >= 9.1

## Keywords
measurement error, straight line fitting, method comparison, allometry