**rags2ridges**
---------------

The R-package **rags2ridges** performs L2-penalized estimation of precison and covariance matrices. The package contains *proper* L2-penalized ML estimators for the precision matrix as well as supporting functions to employ these estimators in a graphical modeling setting. It features fast and efficient algorithms and tools for analysis of relevance networks and gaussian graphical models.

## Installation

The released and tested version of **rags2ridges** is available at
[CRAN](http://cran.r-project.org/package=rags2ridges) (Comprehensive R Archive Network). It can be easily be installed from within R by running

```R
install.packages("rags2ridges")
```

If you wish to install the latest version of **rags2ridges** directly from the master branch here at GitHub, run

```R
#install.packages("devtools")  # Uncomment if devtools is not installed
devtools::install_github("AEBilgrau/rags2ridges")
```

Note, that this version is in development and is different from the version at CRAN. As such, it may be unstable. Be sure that you have the
[package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.

When installed, run `news(package = "rags2ridges")` to view the latest notable changes of GMCM.

For previous versions of **rags2ridges**, visit the old [releases at GitHub](https://github.com/AEBilgrau/rags2ridges/releases) or the [archive at CRAN.](http://cran.r-project.org/src/contrib/Archive/rags2ridges/)


## References

Relevant publications to **rags2ridges** include:

 1. Peeters, C.F.W. and van Wieringen, W.N. (2014) *"rags2ridges: Ridge 
    Estimation of Precision Matrices from High-Dimensional Data. R package"*, 
    version 1.5.
 2. Wessel N. van Wieringen & Carel F.W. Peeters (2014)
    *"Ridge Estimation of Inverse Covariance Matrices from High-Dimensional
    Data"*, arXiv:1403.0904 [stat.ME].
 3. Bilgrau, AE; Peeters CFW; Eriksen, PS; Boegsted, M; & van Wieringen, WN 
    (in preparation) *"Fused Ridge Estimation of Multiple Inverse Covariance 
    Matrices from High-Dimensional Data Classes"*
 4. Carel F.W. Peeters & Wessel N. van Wieringen (in preparation)
    *"The Spectral Condition Number Plot for Regularization Parameter
    Selection"*

Please cite the relevant publication if you use **rags2ridges**.

---
