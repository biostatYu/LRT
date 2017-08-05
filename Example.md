# Example of LRT

## Input data
**[LRT](https://github.com/biostatpzeng/LRT/blob/master/LRT.R)** requires ***y*** (continuous phenotypes, given in the sixth column of [phenotype.fam](https://github.com/biostatpzeng/LRT/blob/master/phenotype.fam)), ***X*** (covariates; including constant 1 for intercept; if no other covariates are present, ***X*** should be given as a n by 1 matrix of 1) and ***G*** (genotypes). It first computes the ***K*** matrix: ***K*** = GG^T (a n by n matrix measureed the similarity of individuals). The rest test is based ON ***y***, ***K*** and ***X***. For the input data, no missing data is allowed. So, missing data should be removed before data analysis.

## [The function of LRT](https://github.com/biostatpzeng/LRT/blob/master/LRT.R)

LRT  <- function(
	y, # phenotype
	X, # covariates
	K, # K matrix
	nsim, # number of values to simulate
	method = "REML" # methods for estimating the mixed model, REML or ML
  )  {

## Testing for the association of cis-SNPs with gene expression level 

***`
library(RLRsim)
library(compiler)
y <- read.table("phenotype.fam")[,6] 
n  <-  length(y)
X  <-  as.matrix(rep(1,n))
k  <- read.table("snp.sXX.txt")
`***

***` fit = LRT (y, X, k, 10e4, method = "REML") `***

fit

$obsLRT

[1] 0.0734309

$pvalue

[1] 0.33444

$lambda

[1] 0.003821497

"$obsLRT" is the LRT statistic, "$pvalue" is the corresponding p-value and "$lambda" is the ratio of sigam_g2/sigam_e2. 




