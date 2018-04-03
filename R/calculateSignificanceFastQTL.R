# Author: Francois Aguet
# Description: This is a post-processing function for FastQTL. It calculates 
# q-values (Storey) and p-value thresholds for all genes in the permutation results file.

suppressMessages(library(qvalue))
suppressMessages(library(tools))
suppressMessages(library(argparser))

# parse inputs
p <- arg_parser("Annotates FastQTL permutation output and runs qvalue")
p <- add_argument(p, "fastqtlOutput", help="")
p <- add_argument(p, "fdr", type="numeric", help="")
p <- add_argument(p, "outfile", help="")
p <- add_argument(p, "--lambda", type="numeric", help="", default=NULL)
args <- parse_args(p)

cat("Processing FastQTL output (", args$fastqtlOutput, ") with FDR=", args$fdr, "\n", sep="")
fastqtl.df <- read.table(args$fastqtlOutput, header=TRUE, stringsAsFactors=FALSE)
stopifnot("pval_beta" %in% colnames(fastqtl.df))

# remove genes w/o variants
nanrows <- is.na(fastqtl.df[, 'pval_beta'])
fastqtl.df <- fastqtl.df[!nanrows, ]
cat("  * Number of genes tested: ", nrow(fastqtl.df), " (excluding ",
    sum(nanrows), " genes w/o variants)\n", sep="")
cat("  * Correlation between Beta-approximated and empirical p-values: ",
    round(cor(fastqtl.df[, 'pval_perm'], fastqtl.df[, 'pval_beta']), 4), "\n", sep="")

# calculate q-values
if (is.null(args$lambda) || is.na(args$lambda)) {
    Q <- qvalue(fastqtl.df[, 'pval_beta'])
} else {
    cat("  * Calculating q-values with lambda = ", args$lambda, "\n", sep="")
    Q <- qvalue(fastqtl.df[, 'pval_beta'], lambda=args$lambda)
}

fastqtl.df$qval <- signif(Q$qvalues, 6)
cat("  * Proportion of significant phenotypes (1-pi0): " , round((1 - Q$pi0), 2), "\n", sep="")
cat("  * eGenes @ FDR ", args$fdr, ":   ", sum(fastqtl.df[, 'qval']<args$fdr), "\n", sep="")

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
ub <- sort(fastqtl.df[fastqtl.df$qval > args$fdr, 'pval_beta'])[1]  # smallest p-value above FDR
lb <- -sort(-fastqtl.df[fastqtl.df$qval <= args$fdr, 'pval_beta'])[1]  # largest p-value below FDR
pthreshold <- (lb+ub)/2
cat("  * min p-value threshold @ FDR ", args$fdr, ": ", pthreshold, "\n", sep="")
fastqtl.df[, 'pval_nominal_threshold'] <- signif(qbeta(pthreshold,
    fastqtl.df[, 'beta_shape1'], fastqtl.df[, 'beta_shape2'], ncp=0, lower.tail=TRUE, log.p=FALSE), 6)

write.table(fastqtl.df, gzfile(args$outfile), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
