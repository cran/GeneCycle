#######################################################################
# Note that this note can directly be run in R.
# Version: GeneCycle 1.0.0 (August 2006)
#######################################################################

#
# EXAMPLE SESSION FOR CELL CYCLE ANALYSIS
# 

# for details see:
#
# Wichert, S., K. Fokianos, and K. Strimmer. 2004.
# Identifying periodically expressed transcripts in microarray
# time series data. Bioinformatics 20:4-20 



# load GeneCycle library
library("GeneCycle")

#######################################################################

# THE DATA:

# the normalized data need to be ready in time series format, i.e. in
# a matrix where each *column* corresponds to a gene, and where the
# *rows* correspond to the individual measurements (time points).

# Example: the Caulobacter data set
data(caulobacter)

# how many samples (11) and how many genes (1444)?
dim(caulobacter)
summary(caulobacter)
get.time.repeats(caulobacter)

# plot first nine time series
plot(caulobacter, 1:9)



#######################################################################

# IDENTIFYING PERIODICALLY EXPRESSED GENES:

# A statistical test developed by Fisher is used to detect 
# periodically expressed genes, and the average periodogram
# is used to visualize the dominant frequencies

# compute and plot average periodogram
avgp.caulobacter <- avgp(caulobacter, "Caulobacter")
avgp.caulobacter

# p-values from Fisher's g test
pval.caulobacter <- fisher.g.test(caulobacter)
pval.caulobacter


#######################################################################
# multiple testing with tail-area based FDR

# test with FDR controlled at on the level 0.05
fdr.out <- fdr.control(pval.caulobacter, Q = 0.05)
fdr.out
fdr.out$num.significant

# proportion of null p-values for different methods
fdr.estimate.eta0(pval.caulobacter, method="conservative")
fdr.estimate.eta0(pval.caulobacter, method="adaptive")
fdr.estimate.eta0(pval.caulobacter, method="bootstrap")
fdr.estimate.eta0(pval.caulobacter, method="smoother")

# FDR test with eta0 != 1
fdr.control(pval.caulobacter, Q = 0.05, eta0=0.9)$num.significant


#######################################################################
# multiple testing with local fdr


# transform exact p-values to z-scores
z.caulobacter <- qnorm(pval.caulobacter)
plot(density.pr(z.caulobacter))

lfdr <- locfdr(z.caulobacter)$fdr
sum(lfdr < 0.2)


#######################################################################


