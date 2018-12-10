###fixing scde based on flemex dep

# in R
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

###reinstall SCDE
require(devtools)
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
