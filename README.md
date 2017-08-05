# LRT: likelihood ratio test (LRT) for SNP set association based on mixed models

**LRT** is a procedure for examining the association of a set of genetic variants (e.g. rare or common SNPs) under the framework of linear mixed models. Very similar to the popular SKAT method which conducts the score test, LRT uses likelihood ratio test to test for the association by testing the variance component parameter. 

Specially, let y be a n by 1 vector of continuous phenotypes on n individuals, X is a n by p covariates, and G is a n by m genotypes of SNPs for a genetic region (i.e., gene). We relate y, X and G by a linear mixed model:

***` y = Xa + Gb + e; b ~ N(0, sigam_g2), e ~　N(0, sigam_ｅ2)　`***

LRT examines the association of G with y (while controling for X) by testing for:

***` H0: b = 0; <==>, H0: sigam_g2 = 0 `***

