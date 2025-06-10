# SPARE

The SPARE package enables recurrent event analysis for GWAS. SPARE has low computational demends, and applies a saddlepoint approximation to achieve statistical accuracy for low-frequency variants. 

The software can be downloaded in R by:
```
library(devtools)
install_github('JasperHof/SPARE')
library(SPARE)
```

The R package enables a GWAS for both .bed and .bgen files using the functions 'SPARE.bed' and 'SPARE.bgen'. Information and examples for running the code can be found by running '?SPARE.bed' and '?SPARE.bgen'. 

A manual for using the SPARE package is included in the repository.

Please do not hesitate to contact me (jasper.hof@qgg.au.dk) or report and issue for any questions or problems in analysis.
