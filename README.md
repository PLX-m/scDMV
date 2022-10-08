# `scDMV`

> A statistical method to detect DNA methylation differences for single-cell bisulfite sequencing data.

__Authors:__ Yan Zhou

## Description

`scDMV` is a statistical method applied to single-cell bisulfite sequencing data(scBS-seq data) to detect differentially methylated regions of DNA. This method is based on a 0-1 inflated  beta binomial distribution model, using the Wald test to calculate p-values for each region in scBS-seq data to identify differentially methylated regions. It outperforms other state-of-the-art methylation difference recognition methods.

## Usage

This GitHub repository is an R package called scDMV, which contains the source code needed to run the `scDMV` method. You need to download ZIP and install the package in R. After downloading the scDMV package, you can directly call the "run_scDMV" function to detect differentially methylated regions. The returns of "run_scDMV" function are p-values and regional level methylation difference (Δ) between two types of samples. In general, regions that satisfy p-values not greater than the p-value cutoff and ∆ greater than the ∆ cutoff are the final DMRs identified by the model. You can specify thresholds to screen for differentially methylated regions based on the p-values and Δ returned by the model.

Here is a description of some of the important parameters in the `run_scDMV` function：

* `treadn`: total reads matrix, the rows represent CpG sites and the columns represent cells.
* `treadx`: methylation reads matrix, the rows represent CpG sites and the columns represent cells.
* `testRegion`: a large list generated by treadn and treadx after CpG sites filtering and region division. Each element of testRegion is also a list, containing the total reads matrix and methylation reads matrix of a region. The rows of the total reads matrix represent CpG sites and the columns represent cells, and so does the methylation reads matrix.
* `sample4c`: the location of the first type sample(e.g 4cell sample), such as sample4c = c(1,2,3,4,5,6,7,8,9,10).
* `sample8c`: the location of the second type sample(e.g 8cell sample), such as sample8c = c(11,12,13,14,15,16).
* `m`: The column where the first gene is, which can be 1.
* `n`: The column where the last gene is, which can be 26.

## Example
```
library(scDMV)
# 8c vs 8c 无差异实验
load("chr3_8cVS8c_data.Rdata")

treadn = chr3_8c8c_data$treadn
treadx = chr3_8c8c_data$treadx
testRegion = chr3_8c8c_data$testRegion
sample8c = c(1:20)
sample4c = c(21:40)

scDMV_result = run_scDMV(treadn,treadx,testRegion,sample8c,sample4c,1,40)
pvalues = scDMV_result$scDMVpvaluem
dlt = scDMV_result$absd

```
'chr3_8cVS8c_data.Rdata' is in the data folder.

## Issues

If you encounter any bugs or have any specific feature requests, please [file an
issue](https://github.com/PLX-m/scDMV/issues).

---


