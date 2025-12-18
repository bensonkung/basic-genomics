# 📊 Summary of Important QC Steps

## Base Data (GWAS) QC:

- Remove missing data - ensures complete SNP information
- Filter ambiguous SNPs (A/T, G/C) - prevents strand alignment issues
- MAF filtering (>0.01) - removes rare variants with unstable estimates
- INFO score filtering (>0.8) - quality control for imputed variants
- Duplicate removal - prevents double-counting of variants
- Effect size transformation - converts OR to log(OR) when needed

## Target Data QC:

- Genotype missingness (<5%) - removes poor quality SNPs
- Sample missingness (<5%) - removes poor quality individuals
- MAF filtering (>0.01) - ensures sufficient allele counts
- HWE test (p>1×10⁻⁶) - detects genotyping errors
- Sex check - identifies sample mix-ups
- Relatedness check - removes related individuals

## LD Clumping:

- Removes correlated variants (r² > 0.1 within 250kb)
- Keeps most significant SNP per LD block
- Critical to avoid double-counting genetic signals

## P-value Thresholds:

- Tests multiple thresholds (5×10⁻⁸ to 1.0)
- Finds "best-fit" that maximizes variance explained
- Balances signal vs. noise
