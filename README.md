## FastQTL

This repository contains a modified version of the FastQTL QTL mapping software, with the following enhancements:

1. Options for filtering by minor allele frequency and minor allele sample count
2. Python wrapper for multi-threaded execution
3. Calculation of q-values (Storey) for FDR estimation (requires R)
4. Minor allele information reported in output

For documentation and the original version, see [QTLtools](https://qtltools.github.io/qtltools/).

Please note: This version of FastQTL is no longer maintained. [tensorQTL](https://github.com/broadinstitute/tensorqtl) is functionally equivalent and implements multiple additional mapping modes.

#### Running multi-threaded analyses

Nominal pass:
```bash
run_FastQTL_threaded.py ${genotypes}.vcf.gz ${phenotypes}.bed.gz ${prefix} --covariates ${covariates}.txt.gz --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 10
```
Permutation pass:
```bash
run_FastQTL_threaded.py ${genotypes}.vcf.gz ${phenotypes}.bed.gz ${prefix} --covariates ${covariates}.txt.gz --permute 1000 10000 --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 10
```
These minor allele filters result in inclusion of genotypes with minor allele frequency ≥ 0.01 and with at least 10 samples containing the minor allele.
