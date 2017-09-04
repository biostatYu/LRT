# LRT: likelihood ratio test for SNP set association based on mixed models

## Introduction
Likelihood ratio test (**LRT**) is a [**R**](https://cran.r-project.org/) procedure for examining the association of a set of genetic variants (e.g. rare or common SNPs) under the framework of linear mixed models. Very similar to the popular [SKAT](http://www.sciencedirect.com/science/article/pii/S0002929711002229?via%3Dihub) method which conducts the score test, LRT uses likelihood ratio test to test for the association by testing the variance component parameter. 

Specifically, let ***y*** be a n by 1 vector of continuous phenotypes on n individuals, ***X*** is a n by p matrix for covariates, and ***G*** is a n by m matrix for genotypes of SNPs for a genetic region (i.e., gene). We relate ***y***, ***X*** and ***G*** by a linear mixed model:

***` y = Xa + Gb + e; b ~ N(0, sigam_g2), e ~　N(0, sigam_e2)　`***

Above, sigam_g2 is the genetic variance, and sigam_e2 is residual variance. 

LRT examines the association of ***G*** with ***y*** (while controling for ***X***) by testing for:

***` H0: b = 0; <==>, H0: sigam_g2 = 0 `***

The LRT statistic is defined as the difference of log-likelihood value between the null model (i.e. ***y = Xa + e***) and the alternative model (i.e. ***y = Xa + Gb + e***). Due to the limited parameter space of sigam_g2, obtaining the null distribution of the LRT statistic is not trivial. We use the method of spectral decomposition proposed by Crainiceanu & Ruppert [(2004)](http://onlinelibrary.wiley.com/wol1/doi/10.1111/j.1467-9868.2004.00438.x/abstract) to generate the exact finite sample null distribution (implemented via the [**R**](https://cran.r-project.org/) package [RLRsim](https://github.com/fabian-s/RLRsim)).


## Note
The LRT procedure was finished in about 2013. At that time, the author (i.e., [Ping Zeng](https://github.com/biostatpzeng)) was a newer and was outside the door of statistical genetics. Thus, you will find that this procedure was not well designed. I put it here for a beautiful recall of that time.

## A new LRT R function, named [ReLRT.R](https://github.com/biostatpzeng/LRT/blob/master/ReLRT.R), was recently rewritten. [ReLRT.R](https://github.com/biostatpzeng/LRT/blob/master/ReLRT.R) esimates the linear mixed models basde on lme function in R package nlme; it also perfrom the approximate LRT.

## Cite
Ping Zeng, Yang Zhao, Jin Liu, Liya Liu, Liwei Zhang, Ting Wang, Shuiping Huang and Feng Chen. Likelihood Ratio Tests in Rare Variant Detection for Continuous Phenotypes. Annals of Human Genetics, 2014, 78(5): 320-332. [DOI: 10.1111/ahg.12071](http://onlinelibrary.wiley.com/wol1/doi/10.1111/ahg.12071/abstract) 

Ciprian M. Crainiceanu and David Ruppert. Likelihood ratio tests in linear mixed models with one variance component. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 2004, 66(1): 165–185. [DOI: 10.1111/j.1467-9868.2004.00438.x](http://onlinelibrary.wiley.com/wol1/doi/10.1111/j.1467-9868.2004.00438.x/abstract) 

Fabian Scheipl, Sonja Greven and Helmut Küchenhoff. Size and power of tests for a zero random effect variance or polynomial regression in additive and linear mixed models. Computational Statistics & Data Analysis, 2008, 52(7): 3283-3299. [DOI: 10.1016/j.csda.2007.10.022](http://www.sciencedirect.com/science/article/pii/S0167947307004306)

R Core Team. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria, 2013 (https://www.R-project.org/). 

## Contact
I am very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn or pingzeng@umich.edu.


