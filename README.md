# LRT: likelihood ratio test (LRT) for SNP set association based on mixed models

## Introduction
**LRT** is a procedure for examining the association of a set of genetic variants (e.g. rare or common SNPs) under the framework of linear mixed models. Very similar to the popular SKAT method which conducts the score test, LRT uses likelihood ratio test to test for the association by testing the variance component parameter. 

Specially, let ***y*** be a n by 1 vector of continuous phenotypes on n individuals, ***X*** is a n by p matrix for covariates, and ***G*** is a n by m matrix for genotypes of SNPs for a genetic region (i.e., gene). We relate y, X and G by a linear mixed model:

***` y = Xa + Gb + e; b ~ N(0, sigam_g2), e ~　N(0, sigam_e2)　`***

Above, sigam_g2 is the genetic variance, and sigam_e2 is residual variance. 

LRT examines the association of G with y (while controling for X) by testing for:

***` H0: b = 0; <==>, H0: sigam_g2 = 0 `***

The LRT statistic is defined as the difference of log-likelihood value between the null model (i.e. y = Xa + e) and the alternative model (i.e. y = Xa + Gb + e). Due to the limited parameter space of sigam_g2, obtaining the null distribution of the LRT statistic is not trivial. We use the method of spectral decomposition proposed by Crainiceanu & Ruppert [(2004)](http://onlinelibrary.wiley.com/wol1/doi/10.1111/j.1467-9868.2004.00438.x/abstract) to generate the exact finite sample null distribution.

## Note
The LRT procedure was finished in about 2013, at that time, the author [Ping Zeng](https://github.com/biostatpzeng) was just outside the door of statistical genetics. So that, this procedure was not well designed. I put it here for a beautiful recall of that time.

## Cite
Ping Zeng, Yang Zhao, Jin Liu, Liya Liu, Liwei Zhang, Ting Wang, Shuiping Huang, Feng Chen. Likelihood Ratio Tests in Rare Variant Detection for Continuous Phenotypes. Annals of Human Genetics, 2014, 78(5): 320-332. 

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn or pingzeng@umich.edu.


