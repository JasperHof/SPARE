# LRM

In LRM, a linear regression on martingale residuals is applied to carry out a recurrent event analysis in GWAS. LRM applies a saddle-point approximation to achieve statistical accuracy for low-frequency variants. 

The software can be downloaded in R by:
```
library(devtools)
install_github('JasperHof/LRM')
library(LRM)
```

The R package enables a GWAS for .bed and .bgen files using the functions 'LRM.bed' and 'LRM.bgen'. Information and examples for running the code can be found by running '?LRM.bed' and '?LRM.bgen'. 

A manual for using the LRM package is included in the repository.

Please do not hesitate to contact me (Jasper.Hof@radboudumc.nl) in case you have any questions or if you encounter a problem in your analysis :-).
