# modSaRa
The modified Screening and Ranking algorithm (modSaRa) can detect chromosome copy number variants with high sensitivity and specificity. For a sequence of intensity values, the modified SaRa will process it by quantile normalization, search for change-point candidates, eliminate unlikely change-points, and then output the potential CNV segments by presenting the start point and end point by SNP or CNV marker index
## Getting Started
The source package needs to be compiled first with the Rtools on Windows. 
## Installing
```
install.packages("devtools")
library(devtools)
install_github("FeifeiXiaoUSC/modSaRa",subdir="package")
```
## Reference Manual
*[Reference Manual](https://github.com/FeifeiXiaoUSC/modSaRa/blob/master/Reference%20manual.pdf)
