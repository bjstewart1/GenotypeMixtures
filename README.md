# GenotypeMixtures
This is an R package which allows users to take the genotype calls from souporcell and stich together genotype clusters across large single cell genomics projects with multiple microfluidic channel experiments where cells are multiplexed in an overlapping mixture design. 

The package aims to simplify the post-processing of experiments initially processed with the souporcell package https://github.com/wheaton5/souporcell (Heaton _et al_. 2020) which allows reference free clustering of cells according to genotype.

### Installation
To install this package use devtools
```
library(devtools)
install_github('bjstewart1/GenotypeMixtures')
```

### Dependencies
* igraph (https://igraph.org/)
* reshape2 (https://cran.r-project.org/web/packages/reshape2/index.html)
* ggplot2 (https://ggplot2.tidyverse.org/)
* pbapply (https://cran.rstudio.com/web/packages/pbapply/index.html)
* ggraph (https://github.com/thomasp85/ggraph)
* vcfR (https://knausb.github.io/vcfR_documentation/)
