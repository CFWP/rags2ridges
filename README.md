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

 1. Peeters, C.F.W., Bilgrau, A.E., & van Wieringen, W.N. (2015). 
    *"rags2ridges: Ridge Estimation of Precision Matrices from High-Dimensional Data"*. 
    R package, version 2.0. 
    URL: https://cran.r-project.org/web/packages/rags2ridges/index.html.
 2. van Wieringen, W.N. & Peeters, C.F.W. (2015).
    *"Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data"*, 
    arXiv:1403.0904v3 [stat.ME].
 3. Bilgrau\*, A.E., Peeters\*, C.F.W., Eriksen, P.S., Boegsted, M., & van Wieringen, W.N. (2015).
    *"Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes"*,
    arXiv:1509.07982v1 [stat.ME]. 
 4. Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2015)
    *"The Spectral Condition Number Plot for Regularization Parameter Determination"*.

Please cite the relevant publications if you use **rags2ridges**.

