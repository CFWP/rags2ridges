![rags2ridges](https://github.com/CFWP/rags2ridges/blob/master/inst/RAGS.png)


[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN version](http://www.r-pkg.org/badges/version/rags2ridges)](http://cran.r-project.org/package=rags2ridges)
[![Coverage Status](https://img.shields.io/codecov/c/github/CFWP/rags2ridges/master.svg)](https://codecov.io/github/CFWP/rags2ridges?branch=master)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/rags2ridges)](http://cran.r-project.org/web/packages/rags2ridges/index.html)
[![Total CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/rags2ridges)](http://www.r-pkg.org/pkg/rags2ridges)



**rags2ridges**
---------------

The R-package **rags2ridges** performs L2-penalized estimation of precison (and covariance) matrices. 
The package contains *proper* L2-penalized ML estimators for the precision matrix as well as supporting functions to employ these estimators in a (integrative or meta-analytic) graphical modeling setting. 
The package has a modular setup and features fast and efficient algorithms.

## Installation

The released and tested version of **rags2ridges** is available at
[CRAN](http://cran.r-project.org/package=rags2ridges) (Comprehensive R Archive Network). It can be easily be installed from within R by running

```R
install.packages("rags2ridges")
```

If you wish to install the latest version of **rags2ridges** directly from the master branch here at GitHub, run

```R
#install.packages("devtools")  # Uncomment if devtools is not installed
devtools::install_github("CFWP/rags2ridges")
```

Note, that this version is in development and is different from the version at CRAN. As such, it may be unstable. Be sure that you have the
[package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.

When installed, run `news(package = "rags2ridges")` to view the latest notable changes to **rags2ridges**.

For previous versions of **rags2ridges**, visit the old [releases at GitHub](https://github.com/AEBilgrau/rags2ridges/releases) or the [archive at CRAN.](http://cran.r-project.org/src/contrib/Archive/rags2ridges/)


## References

Relevant publications to **rags2ridges** include (ordered according to year):

 1. Peeters, C.F.W., Bilgrau, A.E., & van Wieringen, W.N. (2016). 
    *"rags2ridges: Ridge Estimation of Precision Matrices from High-Dimensional Data"*. 
    R package, version 2.1. 
    URL: https://cran.r-project.org/web/packages/rags2ridges/index.html.
 2. van Wieringen, W.N. & Peeters, C.F.W. (To appear).
    *"Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data"*, 
    Computational Statistics & Data Analysis. 
    Available as [arXiv:1403.0904v3 \[stat.ME\]](http://arxiv.org/abs/1403.0904).
 3. Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2016)
    *"The Spectral Condition Number Plot for Regularization Parameter Determination"*.
 4. van Wieringen, W.N. & Peeters, C.F.W. (2015).
    *"Application of a New Ridge Estimator of the Inverse Covariance Matrix
    to the Reconstruction of Gene-Gene Interaction Networks"*.
    In: di Serio, C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.)
    `Computational Intelligence Methods for Bioinformatics and Biostatistics'.
    Lecture Notes in Computer Science, vol. 8623. Springer, pp. 170-179.
 5. Bilgrau\*, A.E., Peeters\*, C.F.W., Eriksen, P.S., Boegsted, M., & van Wieringen, W.N. (2015).
    *"Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes"*.
    Available as [arXiv:1509.07982v1 \[stat.ME\]](http://arxiv.org/abs/1509.07982). 

Please cite the relevant publications if you use **rags2ridges**.
