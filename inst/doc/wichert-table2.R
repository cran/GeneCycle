# last modified: 18 February 2007


# This reproduces Table 2 of Wichert et al. (2004)


# Note: in the following multiple testing is performed in fourways:
#      - as in Wichert et al. using the original Benjamini-Hochberg (1995) algorithm
#         (this yields the "conservative" estimate)
#      - using the FDR approach of Storey and Tibshirani (2003)
#         (this yields the "less conservative" estimate)
#      - using semiparametric Fdr and local fdr, as estimated by the
#         algorithm in fdrtool
#
# Other FDR controlling approaches may lead to different sets of
# significant genes.


library("GeneCycle")
library("qvalue")


# load data sets

data(caulobacter)

# download these files from http://www.strimmerlab.org/data.html 
load("spellmann-yeast.rda.gz")
load("fibroblasts.rda.gz")
load("humanhela.rda.gz")



# calculate p-values and determine number of significant genes
find.periodic.genes <- function(dataset)
{
  cat("Computing p-values ...\n")
  pval = fisher.g.test(dataset)
   
  n1 = sum(qvalue(pval, lambda=0)$qvalues < 0.05)
  n2 = sum(qvalue(pval)$qvalues < 0.05)

  fdr.out <- fdrtool(pval, statistic="pvalue", plot=FALSE)
  n3 <- sum(fdr.out$qval < 0.05) 
  n4 <- sum(fdr.out$lfdr < 0.2) 

  cat("Conservative estimate (Wichert et al.) =", n1, "\n")
  cat("Less conservative estimate =", n2, "\n")
  cat("Semiparametric Fdr < 0.05 (fdrtool) =", n3, "\n")
  cat("Semiparametric fdr < 0.2 (fdrtool) =", n4, "\n")
}


### data analysis ###

dim(cdc15)       # 24 x 4289
find.periodic.genes(cdc15)
# Conservative estimate (Wichert et al.) = 767
# Less conservative estimate = 1293
# Semiparametric Fdr < 0.05 (fdrtool) = 1298 
# Semiparametric fdr < 0.2 (fdrtool) = 1609 

dim(cdc28)       # 17 x 1365
find.periodic.genes(cdc28)
# Conservative estimate (Wichert et al.) = 27   (!! note: typographic error in ms. !!)
# Less conservative estimate = 56
# Semiparametric Fdr < 0.05 (fdrtool) = 67 
# Semiparametric fdr < 0.2 (fdrtool) = 179 

dim(alpha)       # 18 x 4415
find.periodic.genes(alpha)
# Conservative estimate (Wichert et al.) = 469
# Less conservative estimate = 682
# Semiparametric Fdr < 0.05 (fdrtool) = 682 
# Semiparametric fdr < 0.2 (fdrtool) = 976 

dim(elution)     # 14 x 5695
find.periodic.genes(elution)
# Conservative estimate (Wichert et al.) = 194
# Less conservative estimate = 347
# Semiparametric Fdr < 0.05 (fdrtool) = 347 
# Semiparametric fdr < 0.2 (fdrtool) = 764 

dim(caulobacter)     # 11 x 1444
find.periodic.genes(caulobacter)
# Conservative estimate (Wichert et al.) = 45
# Less conservative estimate = 45
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 0 

dim(humanN2)     # 13 x 4574
find.periodic.genes(humanN2)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 1 

dim(humanN3)     # 12 x 5079
find.periodic.genes(humanN3)
# Conservative estimate (Wichert et al.) = 1
# Less conservative estimate = 1
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 0 

dim(score1)      # 12 x 14728
find.periodic.genes(score1)
# Conservative estimate (Wichert et al.) = 1
# Less conservative estimate = 2
# Semiparametric Fdr < 0.05 (fdrtool) = 3 
# Semiparametric fdr < 0.2 (fdrtool) = 37 

dim(score2)      # 26 x 15472
find.periodic.genes(score2)
# Conservative estimate (Wichert et al.) = 135
# Less conservative estimate = 273
# Semiparametric Fdr < 0.05 (fdrtool) = 279 
# Semiparametric fdr < 0.2 (fdrtool) = 978 

dim(score3)      # 48 x 39724
find.periodic.genes(score3)
# Conservative estimate (Wichert et al.) = 6044
# Less conservative estimate = 6971
# Semiparametric Fdr < 0.05 (fdrtool) = 6974 
# Semiparametric fdr < 0.2 (fdrtool) = 7359 

dim(score4)      # 19 x 39192
find.periodic.genes(score4)
# Conservative estimate (Wichert et al.) = 57
# Less conservative estimate = 60
# Semiparametric Fdr < 0.05 (fdrtool) = 60 
# Semiparametric fdr < 0.2 (fdrtool) = 122 

dim(score5)      #  9 x 34890
find.periodic.genes(score5)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0
# Semiparametric Fdr < 0.05 (fdrtool) = 0 
# Semiparametric fdr < 0.2 (fdrtool) = 0 

