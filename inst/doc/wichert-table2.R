
# This reproduces Table 2 of Wichert et al. (2004)


# Note: in the following multiple testing is performed in two ways:
#      - as in Wichert et al. using the original Benjamini-Hochberg (1995) algorithm
#         (this yields the "conservative" estimate)
#      - using the FDR approach of Storey and Tibshirani (2003)
#         (this yields the "less conservative" estimate)
#
# Other FDR controlling approaches may lead to different sets of
# significant genes.


library("GeneCycle")
data(caulobacter)

# download these files from http://www.strimmerlab.org/data.html 
load("spellmann-yeast.rda")
load("fibroblasts.rda")
load("humanhela.rda")

# calculate p-values and determine number of significant genes
find.periodic.genes <- function(dataset)
{
  cat("Computing p-values ...\n")
  pval = fisher.g.test(dataset)
   
  n1 = fdr.control(pval, Q = 0.05, method="conservative", diagnostic.plot=FALSE)$num.significant 
  n2 = fdr.control(pval, Q = 0.05, diagnostic.plot=FALSE)$num.significant

  cat("Conservative estimate (Wichert et al.) =", n1, "\n")
  cat("Less conservative estimate =", n2, "\n")
}


### data analysis ###

dim(cdc15)       # 24 x 4289
find.periodic.genes(cdc15)
# Conservative estimate (Wichert et al.) = 767
# Less conservative estimate = 1293

dim(cdc28)       # 17 x 1365
find.periodic.genes(cdc28)
# Conservative estimate (Wichert et al.) = 27   (!! note: typographic error in ms. !!)
# Less conservative estimate = 56

dim(alpha)       # 18 x 4415
find.periodic.genes(alpha)
# Conservative estimate (Wichert et al.) = 469
# Less conservative estimate = 682

dim(elution)     # 14 x 5695
find.periodic.genes(elution)
# Conservative estimate (Wichert et al.) = 194
# Less conservative estimate = 347

dim(caulobacter)     # 11 x 1444
find.periodic.genes(caulobacter)
# Conservative estimate (Wichert et al.) = 45
# Less conservative estimate = 45

dim(humanN2)     # 13 x 4574
find.periodic.genes(humanN2)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0

dim(humanN3)     # 12 x 5079
find.periodic.genes(humanN3)
# Conservative estimate (Wichert et al.) = 1
# Less conservative estimate = 1

dim(score1)      # 12 x 14728
find.periodic.genes(score1)
# Conservative estimate (Wichert et al.) = 1
# Less conservative estimate = 2

dim(score2)      # 26 x 15472
find.periodic.genes(score2)
# Conservative estimate (Wichert et al.) = 135
# Less conservative estimate = 273

dim(score3)      # 48 x 39724
find.periodic.genes(score3)
# Conservative estimate (Wichert et al.) = 6044
# Less conservative estimate = 6971

dim(score4)      # 19 x 39192
find.periodic.genes(score4)
# Conservative estimate (Wichert et al.) = 57
# Less conservative estimate = 60

dim(score5)      #  9 x 34890
find.periodic.genes(score5)
# Conservative estimate (Wichert et al.) = 0
# Less conservative estimate = 0

