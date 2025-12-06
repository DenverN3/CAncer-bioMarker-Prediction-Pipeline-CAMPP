# CAMPP UK Biobank Examples

This directory contains example templates and workflows for using CAMPP with UK Biobank data.

## Contents

### 1. Example Data Templates

#### `ukb_crc_proteomics_template.txt`
Template for proteomics expression data with UK Biobank formatting:
- **Format:** Tab-delimited text file
- **Structure:**
  - First column: Protein IDs (gene symbols or Uniprot IDs)
  - Subsequent columns: Sample IDs (e.g., UKB_1001, UKB_1002, ...)
  - Values: Normalized Protein Expression (NPX) units from Olink
- **Size:** 20 proteins × 12 samples (template only - expand for real data)

**Note:** This is a minimal template. Actual Olink panels contain 80-96 proteins per panel.

#### `ukb_crc_metadata_template.txt`
Template for sample metadata with all required columns for comprehensive CAMPP analysis:

| Column | Required | Description | Example Values |
|--------|----------|-------------|----------------|
| `ids` | Yes | Sample identifiers (must match proteomics columns) | UKB_1001, UKB_1002 |
| `group` | Yes | Disease status or comparison groups | CRC, Control |
| `batch` | No* | Technical batch identifiers | Batch_A, Batch_B |
| `survival` | No** | Has survival data (1=yes, 0=no for controls) | 0, 1 |
| `age` | No** | Age at baseline (years) | 52, 63 |
| `outcome` | No** | Event status (0=censored, 1=death) | 0, 1 |
| `outcome.time` | No** | Follow-up time (months/years) | 48, 72 |
| `sex` | No*** | Sex (0=Female, 1=Male) | 0, 1 |
| `bmi` | No*** | Body Mass Index | 24.5, 28.3 |
| `smoking` | No*** | Smoking status (0=never, 1=former, 2=current) | 0, 1, 2 |
| `stage` | No*** | Cancer stage (for cases only, NA for controls) | I, II, III, IV, NA |

**Legend:**
- \* Required if using `-b` flag for batch correction
- \** Required if using `-u` flag for survival analysis
- \*** Optional covariates for use with `-r` flag

### 2. Workflow Scripts

#### `run_ukb_crc_analysis.sh`
Complete automated analysis pipeline for UK Biobank colorectal cancer data. This bash script runs all six major CAMPP analysis steps sequentially:

1. **Quality Control & Exploration**
   - Cullen-Frey distribution plots
   - MDS visualization
   - Batch effect correction

2. **Differential Expression Analysis**
   - LIMMA-based statistical testing
   - Covariate adjustment
   - Heatmap generation

3. **LASSO Feature Selection**
   - Cross-validated regularized regression
   - Consensus biomarker identification
   - Classification performance (AUC)

4. **Survival Analysis**
   - Cox proportional hazards models
   - Hazard ratios with 95% CI
   - Forest plot visualization

5. **WGCNA Network Analysis**
   - Co-expression module detection
   - Hub protein identification
   - Network heatmaps

6. **Protein-Protein Interactions**
   - STRING database integration
   - Network arc diagrams

**Usage:**
```bash
# 1. Prepare your data files (must be named exactly):
#    - ukb_crc_proteomics.txt
#    - ukb_crc_metadata.txt

# 2. Run the complete pipeline:
bash examples/run_ukb_crc_analysis.sh

# 3. Results will be saved in Results/ directory
```

**Prerequisites:**
- R version 4.0.0 or higher
- Internet connection (for STRING database access)
- ~5-30 minutes runtime (depends on sample size)

## Quick Start Guide

### Option 1: Test with Template Data

```bash
# Copy templates to main directory
cp examples/ukb_crc_proteomics_template.txt ukb_crc_proteomics.txt
cp examples/ukb_crc_metadata_template.txt ukb_crc_metadata.txt

# Run the analysis pipeline
bash examples/run_ukb_crc_analysis.sh
```

**Note:** Template data is minimal (12 samples) and for demonstration only. Real analyses require adequate sample sizes (≥50 per group recommended).

### Option 2: Use Your Own UK Biobank Data

1. **Obtain UK Biobank proteomics data:**
   - Apply for data access via UK Biobank AMS
   - Download Olink NPX values (Field IDs: 30900-30999)
   - Extract cancer registry linkage data

2. **Format your data:**
   - **Proteomics:** Follow structure of `ukb_crc_proteomics_template.txt`
     - Rows = proteins (gene symbols)
     - Columns = samples (UK Biobank participant IDs)
     - Values = NPX expression levels

   - **Metadata:** Follow structure of `ukb_crc_metadata_template.txt`
     - One row per sample
     - Include all required columns (see table above)
     - Ensure `ids` match proteomics column names exactly

3. **Save files as:**
   - `ukb_crc_proteomics.txt`
   - `ukb_crc_metadata.txt`

4. **Run analysis:**
   ```bash
   bash examples/run_ukb_crc_analysis.sh
   ```

### Option 3: Custom Manual Analysis

For fine-grained control, run individual CAMPP steps:

```bash
# Example: Run only differential expression
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n MyAnalysis \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -r TRUE,age,sex,bmi \
  -a DE
```

See `UK_BIOBANK_USE_CASE.md` for detailed parameter explanations.

## Customization Tips

### Adjusting Significance Cutoffs

Edit the `-f` parameter in the script:
```bash
-f 0.5,0.05   # logFC > 0.5, FDR < 0.05 (default)
-f 1.0,0.01   # More stringent: logFC > 1.0, FDR < 0.01
-f 0.3,0.10   # More lenient: logFC > 0.3, FDR < 0.10
```

### Changing LASSO to Elastic Net

Edit the `-l` parameter:
```bash
-l 1.0    # LASSO (alpha=1.0)
-l 0.9    # Elastic Net with 90% L1 penalty
-l 0.5    # Elastic Net with equal L1/L2 penalty
```

### Adding/Removing Covariates

Edit the `-r` parameter in Steps 2 and 4:
```bash
# Include covariates in both DE and survival:
-r TRUE,age,sex,bmi,smoking

# Include covariates ONLY in survival analysis:
-r FALSE,age,sex,bmi,smoking

# Different covariates for DE vs survival:
# Step 2 (DE): -r TRUE,age,sex
# Step 4 (Survival): -r TRUE,age,sex,stage,treatment
```

### Modifying WGCNA Parameters

Edit the `-x` parameter:
```bash
-x 5,25,30    # Min module size=5, merge at 25% dissimilarity, top 30% hubs
-x 10,20,50   # Larger modules, more stringent merging, more hubs returned
-x 3,30,25    # Smaller modules, lenient merging, fewer hubs
```

### Changing Color Scheme

Edit the `-c` parameter (must match number of groups):
```bash
-c darkred,steelblue     # 2 groups (CRC vs Control)
-c red,blue,green        # 3 groups
-c "#FF5733","#33FF57"   # Hex color codes
```

See R color reference: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

## Expected Runtime

| Analysis Step | Sample Size | Approximate Time |
|---------------|-------------|------------------|
| QC + MDS | 300 samples | 2-3 minutes |
| Differential Expression | 300 samples | 3-5 minutes |
| LASSO | 300 samples, 92 proteins | 5-7 minutes |
| Survival Analysis | 150 cases | 2-4 minutes |
| WGCNA | 50 proteins | 5-10 minutes |
| PPI Network | 20-30 DE proteins | 3-5 minutes |
| **Total** | **~20-35 minutes** |

**Note:** Runtime increases with larger sample sizes and feature counts. WGCNA is particularly sensitive to feature count (not recommended for >5000 features).

## Troubleshooting

### Common Issues

**1. "Error: could not find function"**
- **Cause:** Missing R package dependencies
- **Solution:** Let CAMPP auto-install packages on first run, or use `-e stable` flag for renv library

**2. "Error: ids in metadata do not match data columns"**
- **Cause:** Mismatch between sample IDs in proteomics and metadata
- **Solution:** Ensure `ids` column in metadata exactly matches column names in proteomics file

**3. "Warning: batch is confounded with group"**
- **Cause:** All cases in some batches, all controls in others
- **Solution:** Batch correction may not work properly. Redesign study with balanced batches.

**4. "Error: too few samples for LASSO"**
- **Cause:** <9 samples in one of the groups
- **Solution:** Increase sample size or skip LASSO step

**5. "STRING database connection failed"**
- **Cause:** No internet connection or server down
- **Solution:** Check connectivity, or skip PPI analysis (`-p` flag)

### Getting Help

- **Full documentation:** See `../UK_BIOBANK_USE_CASE.md`
- **CAMPP Manual:** See `../CAMPPManual.pdf`
- **CAMPP help:** Run `Rscript CAMPP.R -h`
- **Contact authors:** thilde.terkelsen@sund.ku.dk, elenap@cancer.dk

## Data Privacy and Ethics

**IMPORTANT:** UK Biobank data is subject to strict data protection requirements.

### Data Access Requirements
1. Approved UK Biobank research application
2. Institutional ethics approval
3. Data Transfer Agreement in place
4. Secure data storage (encrypted, access-controlled)

### What NOT to Share
- **Never commit actual UK Biobank data to GitHub**
- **Never share participant IDs publicly**
- **Never publish data that could re-identify individuals**

### What Can Be Shared
✓ Template data structures (like those in this directory)
✓ Analysis code and workflows
✓ Aggregated results (summary statistics, plots)
✓ Documentation and use cases

For more information:
- UK Biobank Access Management System: https://www.ukbiobank.ac.uk/
- UK Biobank Ethics: https://www.ukbiobank.ac.uk/learn-more-about-uk-biobank/about-us/ethics

## Citation

If you use these templates and workflows in your research, please cite:

**CAMPP:**
Terkelsen, T., Krogh, A., & Papaleo, E. (2020). CAncer bioMarker Prediction Pipeline (CAMPP)—A standardized framework for the analysis of quantitative biological data. *PLoS Computational Biology*, 16(3), e1007665.

**UK Biobank:**
Sudlow, C., et al. (2015). UK biobank: an open access resource for identifying the causes of a wide range of complex diseases of middle and old age. *PLoS Medicine*, 12(3), e1001779.

## License

These examples are provided under the same license as CAMPP. See `../LICENSE.md`.

UK Biobank data use is governed by the UK Biobank Material Transfer Agreement.
