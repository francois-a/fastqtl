# Author: Francois Aguet
# Description: This is a post-processing function for FastQTL. It calculates 
# q-values (Storey) and p-value thresholds for all genes in the permutation results file.

suppressMessages(library(qvalue))
suppressMessages(library(tools))

# parse inputs
args <- commandArgs(trailingOnly=TRUE)
fastqtlOutput <- args[1]
fdr <- as.numeric(args[2])
outfile <- args[3]

cat("Processing FastQTL output (", fastqtlOutput, "), with FDR=", fdr, "\n", sep="")

# input files have no headers
D <- read.table(fastqtlOutput, header=FALSE, stringsAsFactors=FALSE)
if (dim(D)[2]==17) {
    colnames(D) <- c('gene_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'variant_id', 'tss_distance',
        'minor_allele_samples', 'minor_allele_count', 'maf', 'ref_factor',
        'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta')
} else {
    stop("FastQTL output in unrecognized format (mismatched number of columns).")
}

# remove genes w/o variants
nanrows <- is.na(D[, 'pval_beta'])
D <- D[!nanrows, ]
cat("  * Number of genes tested: ", nrow(D), " (excluding ", sum(nanrows), " genes w/o variants)\n", sep="")
cat("  * Correlation between Beta-approximated and empirical p-values: ", round(cor(D[, 'pval_perm'], D[, 'pval_beta']), 4), "\n", sep="")

# calculate q-values
Q <- qvalue(D[, 'pval_beta'])
D$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", fdr, ":   ", sum(D[, 'qval']<0.05), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(D[D$qval > fdr, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-D[D$qval <= fdr, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("  * min p-value threshold @ FDR ", fdr, ": ", pthreshold, "\n", sep="")
D[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold, D[, 'beta_shape1'], D[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

write.table(D, gzfile(outfile), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
