
# rags2ridges pkgdown-website

This is a small README for the **rags2ridges** documentation website created by the **pkgdown**-package.
The `_pkgdown.yml` file primarily controls how the site is build. 
Learn more at: https://pkgdown.r-lib.org/


# Local building of site

To build the site run:

```r
install.packages("pkgdown")
pkgdown::build_site()
```

Individual sections of the site can be build in a similar manner.
After successfully building, a `/docs` folder containing the generated website will appear. These are local files only and should be committed.


# Online site

Travis-CI will be responsible for building (generating `/docs`)and pushing it to a specific git(hub) branch. When build, GitHub hosts the site and exposes it to the internet (via 'GitHub pages').
