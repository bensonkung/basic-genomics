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

# Quick Start Guide - PRS Pipeline

## For Standard PRS Analysis

### 1. Edit Configuration
```bash
vim prs_config.sh  # or your other favorite text editor
```

Set these required parameters:
- `TARGET_GENO` = path to your genotype files (without extension)
- `GWAS_SUMSTATS` = path to GWAS summary statistics
- `SNP_COL`, `A1_COL`, `BETA_COL`, `PVAL_COL` = column numbers in GWAS file
- `OUTPUT_DIR` and `PROJECT_NAME`

### 2. Run Pipeline
```bash
./calculate_prs.sh
```

### 3. Analyze Results
```bash
# Replace with your phenotype file and column names
Rscript prs_analysis.R my_prs_analysis phenotype.txt phenotype_column PC1,PC2,PC3
```

## Expected Runtime

| Data Size | Runtime |
|-----------|---------|
| Small (1K samples, 100K SNPs) | 15-30 min |
| Medium (10K samples, 1M SNPs) | 1-2 hours |
| Large (100K samples, 10M SNPs) | 4-8 hours |
| Cell-type (461 clusters) | 2-4 hours |

## Common First-Time Issues

**"Column out of bounds"**
→ Check your SNP_COL, A1_COL, etc. in config
→ Open GWAS file and count columns (starting from 1)

**"No variants found"**
→ Check SNP ID format matches (rs IDs vs chr:pos)
→ Check genome build matches (hg19 vs hg38)

**"Permission denied"**
→ Run: `chmod +x *.sh`

**"plink2: command not found"**
→ On Sherlock, enter "ml biology" then "ml plink/2.0a7"


## Next Steps After First Run

1. **Check the log file** for any warnings
2. **Review QC report** (`*_QC_REPORT.txt`)
3. **Test association** with phenotype
4. **Adjust parameters** if needed
5. **Run validation** in independent sample

---

## Getting Help

1. Check `README.md` for detailed documentation
2. Review log files for specific errors
3. Check PLINK2 docs: https://www.cog-genomics.org/plink/2.0/
4. PRS Tutorial: https://choishingwan.github.io/PRS-Tutorial/
