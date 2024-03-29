
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN version](https://www.r-pkg.org/badges/version/rags2ridges)](https://cran.r-project.org/package=rags2ridges)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/rags2ridges)](https://cran.r-project.org/package=rags2ridges/index.html)
[![Total CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/rags2ridges)](https://www.r-pkg.org/pkg/rags2ridges)

<img src="man/figures/RAGS.png" alt="" />

**rags2ridges**
---------------

The R-package **rags2ridges** performs L2-penalized estimation of precison (and covariance) matrices. 
The package contains *proper* L2-penalized maximum-likelihood estimators for the precision matrix as well as supporting functions to employ these estimators in a (integrative or meta-analytic) graphical modeling setting. The package has a modular setup and features fast and efficient algorithms.


## Installation

The released and tested version of **rags2ridges** is available at
[CRAN](https://cran.r-project.org/package=rags2ridges) (Comprehensive R Archive Network). It can be easily be installed from within R by running

```R
install.packages("rags2ridges")
```

If you wish to install the latest version of **rags2ridges** directly from the master branch here at GitHub, run

```R
#install.packages("remotes")  # Uncomment if not installed
remotes::install_github("CFWP/rags2ridges")
```

Note, that this version is in development and is different from the version at CRAN. As such, it may be unstable. Be sure that you have the
[package development prerequisites](https://support.posit.co/hc/en-us/articles/200486498-Package-Development-Prerequisites) if you wish to install the package from the source.

Visit [CRAN](https://cran.r-project.org/package=rags2ridges/news/news.html), [the documentation site](https://cfwp.github.io/rags2ridges/news/index.html), or run `news(package = "rags2ridges")` after installation to view the latest notable changes to **rags2ridges**. 

For previous versions of **rags2ridges**, visit the [archive at CRAN.](https://cran.r-project.org/src/contrib/Archive/rags2ridges/)


## Usage and getting started

The `vignette("rags2ridges")` provides a light introduction to **rags2ridges** and details how to quickly get started.


## References

Relevant publications to **rags2ridges** include (ordered according to year):

 1. Peeters, C.F.W., Bilgrau, A.E., & van Wieringen, W.N. (2022). 
    *"rags2ridges: A One-Stop-l2-Shop for Graphical Modeling of High-Dimensional Precision Matrices"*. 
    Journal of Statistical Software, vol. 102(4):1-32.
    ([doi:10.18637/jss.v102.i04](https://doi.org/10.18637/jss.v102.i04)).
 2. Peeters, C.F.W., van de Wiel, M.A., & van Wieringen, W.N. (2020)
    *"The Spectral Condition Number Plot for Regularization Parameter Evaluation"*,
    Computational Statistics, vol. 35:629-646
    ([doi:10.1007/s00180-019-00912-z](https://doi.org/10.1007/s00180-019-00912-z)).
 3. Bilgrau\*, A.E., Peeters\*, C.F.W., Eriksen, P.S., Boegsted, M., & van Wieringen, W.N. (2020).
    *"Targeted Fused Ridge Estimation of Inverse Covariance Matrices from Multiple High-Dimensional Data Classes"*,
    Journal of Machine Learning Research, vol. 21(26):1-52
    ([PDF](https://www.jmlr.org/papers/volume21/15-509/15-509.pdf)).
 4. van Wieringen, W.N. & Peeters, C.F.W. (2016).
    *"Ridge Estimation of Inverse Covariance Matrices from High-Dimensional Data"*, 
    Computational Statistics & Data Analysis, vol. 103:284-303
    ([doi:10.1016/j.csda.2016.05.012](https://www.sciencedirect.com/science/article/pii/S0167947316301141)).
 5. van Wieringen, W.N. & Peeters, C.F.W. (2015).
    *"Application of a New Ridge Estimator of the Inverse Covariance Matrix
    to the Reconstruction of Gene-Gene Interaction Networks"*.
    In: di Serio, C., Lio, P., Nonis, A., and Tagliaferri, R. (Eds.)
    `Computational Intelligence Methods for Bioinformatics and Biostatistics'.
    Lecture Notes in Computer Science, vol. 8623. Springer, pp. 170-179
    ([doi:10.1007/978-3-319-24462-4_15](https://link.springer.com/chapter/10.1007/978-3-319-24462-4_15)).

Please cite the relevant publications if you use **rags2ridges**.
