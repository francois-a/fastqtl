# Author: Francois Aguet
# Description: This is a post-processing function for FastQTL. It calculates 
# q-values (Storey) and p-value thresholds for all genes in the permutation results file.

suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(argparser))

# parse inputs
p <- arg_parser("Annotates FastQTL permutation output and runs qvalue")
p <- add_argument(p, "fastqtlOutput", help="")
p <- add_argument(p, "fdr", help="")
p <- add_argument(p, "outfile", help="")
p <- add_argument(p, "--lambda", type="numeric", help="", default=NULL)
args <- parse_args(p)

cat("Processing FastQTL output (", args$fastqtlOutput, "), with FDR=", args$fdr, "\n", sep="")

# input files have no headers
D <- read.table(args$fastqtlOutput, header=FALSE, stringsAsFactors=FALSE)
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
if (is.null(args$lambda) || is.na(args$lambda)) {
    Q <- qvalue(D[, 'pval_beta'])
} else {
    cat("  * Calculating q-values with lambda = ", args$lambda, "\n", sep="")
    Q <- qvalue(D[, 'pval_beta'], lambda=args$lambda)
}

D$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", args$fdr, ":   ", sum(D[, 'qval']<args$fdr), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(D[D$qval > args$fdr, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-D[D$qval <= args$fdr, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("  * min p-value threshold @ FDR ", args$fdr, ": ", pthreshold, "\n", sep="")
D[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold, D[, 'beta_shape1'], D[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

write.table(D, gzfile(args$outfile), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
