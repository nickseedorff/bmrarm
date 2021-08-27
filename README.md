# bmrarm

<!-- badges: start -->
<!-- badges: end -->

`bmrarm` A package for fitting longitudinal and time series models with mixed outcome type (continuous and ordinal).

## Installation

### Prerequisites
```{r eval = FALSE}
install.packages("devtools")
library("devtools")
```

### Without vignette
```{r eval = FALSE}
## Download the package without the vignette
install_github("nickseedorff/bmrarm")
library(bmrarm)
```

### With vignette
```{r eval = FALSE}
install_github("nickseedorff/bmrarm", build_vignettes = T)
library(bmrarm)

## View the vignette
vignette("Introduction", package = "bmrarm")
```
