# UK Biobank Use Case: Colorectal Cancer Biomarker Discovery

## Overview

This use case demonstrates how to use the CAMPP (CAncer bioMarker Prediction Pipeline) with UK Biobank data for discovering and validating serum protein biomarkers for colorectal cancer (CRC) risk stratification.

## Scientific Background

### Research Question
Can we identify serum protein biomarkers from UK Biobank participants that:
1. Differentiate colorectal cancer patients from healthy controls
2. Predict overall survival in CRC patients
3. Show coordinated expression patterns via network analysis

### UK Biobank Context
The UK Biobank is a large-scale prospective cohort study with over 500,000 participants aged 40-69 years recruited between 2006-2010. It includes:
- Baseline health assessments and biological samples
- Long-term health outcome tracking via NHS record linkage
- Proteomics data (Olink panels) for ~50,000 participants
- Cancer registry linkage for incident cancer diagnoses

## Use Case Scenario

### Cohort Selection

**Cases (n=150):**
- UK Biobank participants diagnosed with colorectal cancer (ICD-10: C18-C20)
- Serum samples collected at baseline (pre-diagnosis)
- Complete follow-up data available (5+ years)
- Olink proteomics data available (Oncology II panel, ~90 proteins)

**Controls (n=150):**
- Age-matched (±2 years) healthy participants
- No cancer diagnosis during follow-up period
- Matched for sex, BMI category, and smoking status
- Same proteomics platform and batch processing

**Technical Batches:**
- Samples processed across 6 different Olink assay plates (Batch_A through Batch_F)
- Each batch contains balanced case/control ratios

### Data Structure

#### 1. Proteomics Expression Data (`ukb_crc_proteomics.txt`)
- **Rows:** 92 proteins (Olink Oncology II panel)
- **Columns:** 300 samples (150 cases + 150 controls)
- **Values:** Normalized Protein Expression (NPX) units
- **Format:** Tab-delimited text file

```
IDs         UKB_1001    UKB_1002    UKB_1003    ...
CEACAM5     5.23        6.78        5.12        ...
CA125       3.45        4.21        3.89        ...
AFP         2.11        2.34        2.08        ...
...
```

#### 2. Metadata (`ukb_crc_metadata.txt`)
Required columns for CAMPP analysis:

| Column | Description | Example Values |
|--------|-------------|----------------|
| `ids` | UK Biobank participant ID | UKB_1001, UKB_1002 |
| `group` | Case/control status | CRC, Control |
| `batch` | Olink assay plate | Batch_A, Batch_B, ... Batch_F |
| `survival` | Has survival data (1=yes, 0=no) | 1 (for CRC), 0 (for Control) |
| `age` | Age at baseline (years) | 52, 63, 58 |
| `outcome` | Death status (1=death, 0=censored) | 0, 1 |
| `outcome.time` | Follow-up time (months) | 72, 48, 96 |
| `sex` | Sex (0=Female, 1=Male) | 0, 1 |
| `bmi` | Body Mass Index (kg/m²) | 24.5, 28.3, 31.2 |
| `smoking` | Smoking status (0=never, 1=former, 2=current) | 0, 1, 2 |
| `stage` | Cancer stage (for cases only) | I, II, III, IV, NA |

Example structure:
```
ids         group    batch     survival  age  outcome  outcome.time  sex  bmi   smoking  stage
UKB_1001    CRC      Batch_A   1         62   1        48            1    27.3  1        III
UKB_1002    Control  Batch_A   0         61   0        84            1    26.8  1        NA
UKB_1003    CRC      Batch_B   1         58   0        72            0    24.1  0        II
...
```

## Analysis Workflow

### Step 1: Data Quality Control and Preprocessing

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_Discovery \
  -z median \
  -t log2 \
  -b batch \
  -j TRUE \
  -s TRUE \
  -c darkred,steelblue \
  -k group
```

**Parameters Explained:**
- `-v other`: Proteomics NPX data (not array/seq/ms)
- `-z median`: Median centering for normalization
- `-t log2`: Log2 transformation of NPX values
- `-b batch`: Correct for Olink plate batch effects using ComBat
- `-j TRUE`: Generate Cullen-Frey distribution plots for QC
- `-s TRUE`: Create MDS plot for sample clustering visualization
- `-c darkred,steelblue`: Custom colors for CRC (red) and Control (blue)
- `-k group`: Label MDS plot points by case/control status

**Expected Outputs:**
- `UKB_CRC_Discovery_MDSplot.pdf`: Sample clustering before analysis
- `DataChecks/*.pdf`: Distribution plots for 10 random proteins
- `CAMPPlog.txt`: Processing log with QC metrics

---

### Step 2: Differential Expression Analysis

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_DE \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -r TRUE,age,sex,bmi,smoking \
  -a DE \
  -s TRUE
```

**Parameters Explained:**
- `-f 0.5,0.05`: Significance cutoffs (|logFC| > 0.5, FDR < 0.05)
- `-r TRUE,age,sex,bmi,smoking`: Include covariates in DE model
  - `TRUE` = use covariates in both DE analysis and survival analysis
- `-a DE`: Generate heatmap of differentially expressed proteins

**Expected Outputs:**
- `Results/DEAAResults/UKB_CRC_DE_DE.txt`: DE results with logFC, p-values, FDR
- `UKB_CRC_DE_heatmap.pdf`: Clustered heatmap of significant proteins
- `UKB_CRC_DE_MDSplot.pdf`: Sample visualization

**Interpretation:**
- Identify proteins significantly altered in pre-diagnostic CRC serum
- Control for confounders (age, sex, BMI, smoking)
- Batch-corrected expression ensures technical variability doesn't drive results

---

### Step 3: LASSO Regularized Regression for Feature Selection

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_LASSO \
  -z median \
  -t log2 \
  -b batch \
  -l 1.0 \
  -a LASSO
```

**Parameters Explained:**
- `-l 1.0`: LASSO regression (alpha=1.0)
  - For Elastic Net, use values like 0.5, 0.7, 0.9
- `-a LASSO`: Generate heatmap of LASSO-selected proteins

**Expected Outputs:**
- `Results/LASSOResults/UKB_CRC_LASSO_LASSO.txt`: Selected protein panel
- `Results/LASSOResults/UKB_CRC_LASSO_CrossValidationPlot.pdf`: CV error curves
- `Results/LASSOResults/UKB_CRC_LASSO_AUC.txt`: Classification performance (AUC)
- `Results/LASSOResults/UKB_CRC_LASSO_overlap_DEAA_LASSO_EN.pdf`: Venn diagram
- `Results/LASSOResults/UKB_CRC_LASSO_DEA_LASSO_Consensus.txt`: Consensus biomarkers

**Interpretation:**
- LASSO selects minimal protein set with optimal predictive power
- Cross-validation ensures model stability across random splits
- Consensus proteins appear in both DE and LASSO = high-confidence biomarkers
- AUC quantifies diagnostic potential (AUC > 0.80 = strong classifier)

---

### Step 4: Survival Analysis with Cox Proportional Hazards

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_Survival \
  -z median \
  -t log2 \
  -b batch \
  -u Consensus \
  -r TRUE,age,sex,bmi,smoking,stage \
  -q 20
```

**Parameters Explained:**
- `-u Consensus`: Test consensus proteins from DE ∩ LASSO
  - Alternatives: `DE`, `LASSO`, `EN`, `ALL`
- `-r TRUE,age,sex,bmi,smoking,stage`: Covariates in Cox model
- `-q 20`: Plot max 20 proteins per forest plot page

**Expected Outputs:**
- `Results/SurvivalResults/UKB_CRC_Survival_survival.txt`: Hazard ratios, 95% CI, p-values
- `Results/SurvivalResults/UKB_CRC_Survival_individual_corrplots.pdf`: Forest plots

**Interpretation:**
- Hazard Ratio (HR) > 1: High protein expression → worse survival
- HR < 1: High protein expression → better survival (protective)
- Adjusted for clinical covariates (stage is critical confounder)
- Identifies prognostic biomarkers within CRC cases

**IMPORTANT NOTE:**
The pipeline automatically tests proportional hazards assumptions. If covariates fail the PH test, you may need to stratify:

```bash
# If 'smoking' violates PH assumption:
Rscript CAMPP.R \
  [... same parameters as above ...] \
  -y smoking
```

---

### Step 5: Weighted Gene Co-expression Network Analysis (WGCNA)

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_WGCNA \
  -z median \
  -t log2 \
  -b batch \
  -w DE \
  -x 5,25,30
```

**Parameters Explained:**
- `-w DE`: Use only differentially expressed proteins for WGCNA
  - Alternatives: `DA`, `ALL` (all proteins - not recommended if >5000 features)
- `-x 5,25,30`: WGCNA cutoffs
  1. Minimum module size = 5 proteins
  2. Module merge threshold = 25% dissimilarity
  3. Return top 30% most connected proteins per module

**Expected Outputs:**
- `Results/WGCNAResults/WGCNAPlots/UKB_CRC_WGCNA_softpowerplot.pdf`: Power selection
- `Results/WGCNAResults/WGCNAPlots/UKB_CRC_WGCNA_moduleTree.pdf`: Dendrogram with colors
- `Results/WGCNAResults/WGCNAPlots/UKB_CRC_WGCNA_module[color]_moduleHM.pdf`: Network heatmaps
- `Results/WGCNAResults/WGCNAPlots/UKB_CRC_WGCNA_[color]_moduleIC.pdf`: Hub protein rankings
- `Results/WGCNAResults/module[color]_moduleRes.txt`: Module membership + IC scores

**Interpretation:**
- Identifies co-expressed protein modules
- Hub proteins (high IC) = master regulators or pathway representatives
- Modules may represent distinct biological processes
- Can correlate modules with clinical traits (stage, survival time)

---

### Step 6: Protein-Protein Interaction Networks

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_PPI \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -p hgnc_symbol,11.0
```

**Parameters Explained:**
- `-p hgnc_symbol,11.0`:
  1. Gene identifier type (must match protein names in data)
  2. STRING database version 11.0

**Alternative Identifiers:**
- `ensembl_peptide_id`
- `ensembl_gene_id`
- `ensembl_transcript_id`
- `uniprotswissprot`

**Expected Outputs:**
- `Results/InteractionResults/CRCvsControl_AllInteractions.txt`: Edge list with confidence scores
- `Results/InteractionResults/InteractionPlots/CRCvsControl_TopInteractions.tiff`: Arc diagram

**Interpretation:**
- Maps DE proteins onto STRING PPI network
- High-confidence edges (score > 700) indicate known interactions
- Network hubs may be therapeutic targets
- Identifies functionally connected biomarker clusters

---

### Step 7: K-means Clustering for Sample Subtyping

```bash
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n UKB_CRC_Kmeans \
  -z median \
  -t log2 \
  -b batch \
  -k group
```

**Parameters Explained:**
- `-k group`: Color/label MDS points by case/control status
  - Can use any metadata column (e.g., `stage`, `smoking`)

**Expected Outputs:**
- `Results/KmeansResults/BestKmeans_C[k].pdf`: MDS plot with cluster assignments
- `Results/KmeansResults/UKB_CRC_Kmeans_Metadata_Kmeans.txt`: Cluster membership

**Interpretation:**
- Identifies molecular subtypes within CRC cases
- BIC-optimized cluster number (automatic selection)
- Can test if clusters associate with survival, stage, or response

---

## Comprehensive Multi-Step Workflow

For a complete analysis pipeline combining all steps:

```bash
#!/bin/bash
# UK Biobank CRC Biomarker Discovery - Complete Pipeline

PROJECT="UKB_CRC_Complete"

# Step 1: QC and exploration
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n ${PROJECT}_QC \
  -z median \
  -t log2 \
  -b batch \
  -j TRUE \
  -s TRUE \
  -c darkred,steelblue \
  -k group

# Step 2: Differential expression with covariates
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n ${PROJECT}_DE \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -r TRUE,age,sex,bmi,smoking \
  -a DE

# Step 3: LASSO feature selection
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n ${PROJECT}_LASSO \
  -z median \
  -t log2 \
  -b batch \
  -l 1.0 \
  -a Consensus

# Step 4: Survival analysis on consensus biomarkers
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n ${PROJECT}_Survival \
  -z median \
  -t log2 \
  -b batch \
  -u Consensus \
  -r TRUE,age,sex,bmi,smoking,stage \
  -q 20

# Step 5: WGCNA network analysis
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n ${PROJECT}_WGCNA \
  -z median \
  -t log2 \
  -b batch \
  -w DE \
  -x 5,25,30

# Step 6: PPI networks
Rscript CAMPP.R \
  -d ukb_crc_proteomics.txt \
  -m ukb_crc_metadata.txt \
  -v other \
  -g ids,group \
  -n ${PROJECT}_PPI \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -p hgnc_symbol,11.0

echo "Analysis complete! Results in Results/ directory"
```

---

## Expected Results Summary

### Key Findings from Hypothetical Analysis

**1. Differential Expression:**
- 23 proteins significantly altered in pre-diagnostic CRC serum
- Top hits: CEA, CA19-9, MMP7, VEGF-A, IL-6
- logFC ranges: 0.6 to 2.3 (upregulated in cases)

**2. LASSO Selection:**
- Optimal 8-protein signature for CRC classification
- Cross-validated AUC: 0.84 (95% CI: 0.79-0.89)
- Consensus with DE: 6 proteins

**3. Survival Analysis:**
- High IL-6: HR=2.1, p<0.001 (poor prognosis)
- High TIMP1: HR=1.8, p=0.003 (poor prognosis)
- High OPG: HR=0.6, p=0.012 (protective)

**4. WGCNA:**
- 3 major modules identified:
  - Red module: Inflammation/immune response (IL-6, TNF-α, IL-8)
  - Blue module: ECM remodeling (MMPs, TIMPs)
  - Green module: Angiogenesis (VEGF-A, ANGPT2)

**5. PPI Networks:**
- Dense cluster around IL-6/JAK-STAT pathway
- MMP-TIMP interactions highly enriched
- Novel connection: CEACAM5-MUC16 interaction

---

## Data Availability and Access

### UK Biobank Data Access
1. Register as an approved UK Biobank researcher
2. Submit project application for proteomics data access
3. Download data via UK Biobank data portal
4. Proteomics data fields:
   - Field ID: 30900-30999 (Olink panel NPX values)
   - Cancer registry linkage: Available for all participants

### Example Data Format
This repository includes example template files:
- `examples/ukb_crc_proteomics_template.txt`
- `examples/ukb_crc_metadata_template.txt`

**Note:** These are templates only. Actual UK Biobank data must be obtained through official channels.

---

## Quality Control Recommendations

### Before Running CAMPP

1. **Sample Size Requirements:**
   - Minimum 9 samples per group for LASSO
   - Survival analysis: adequate events (deaths) for Cox models
   - Recommended: ≥50 samples per group for stable estimates

2. **Missing Data:**
   - CAMPP auto-removes features with >70% missing values
   - CAMPP auto-removes samples with >80% missing values
   - Remaining NAs imputed via LLS/KNN

3. **Batch Effects:**
   - Always check for batch confounding with group
   - Use `-s TRUE` to visualize batch effects in MDS
   - ComBat correction (`-b batch`) recommended for multi-plate studies

4. **Outlier Detection:**
   - Review MDS plots for extreme outliers
   - Check CAMPPlog.txt for QC warnings
   - Consider excluding samples >4 SD from group centroid

### After Running CAMPP

1. **Validation:**
   - External cohort replication (essential for biomarker discovery)
   - Bootstrap resampling for stability assessment
   - ROC curve validation on held-out test set

2. **Biological Interpretation:**
   - Pathway enrichment analysis (KEGG, Reactome)
   - Literature review of top biomarkers
   - Correlation with known cancer markers (CEA, CA19-9)

3. **Clinical Utility:**
   - Calculate Net Reclassification Improvement (NRI)
   - Decision curve analysis for clinical benefit
   - Cost-effectiveness assessment

---

## Troubleshooting

### Common Issues

**1. "Insufficient samples for LASSO"**
- Ensure ≥9 samples per group
- Check metadata `group` column for typos

**2. "Proportional hazards assumption violated"**
- Use `-y` flag to stratify problematic covariates
- Example: `-y smoking,stage`

**3. "WGCNA fails to find modules"**
- Try lowering minimum module size: `-x 3,25,30`
- Ensure using DE proteins (`-w DE`) not all features

**4. "STRING database connection error"**
- Check internet connectivity
- Ensure gene identifiers match `-p` specification
- Verify protein names don't have special characters

**5. "ComBat batch correction warning"**
- Ensure batches aren't perfectly confounded with groups
- Check batch distribution across cases/controls

---

## Citation

If using this workflow with UK Biobank data, please cite:

1. **CAMPP Publication:**
   Terkelsen, Thilde, Anders Krogh, and Elena Papaleo. "CAncer bioMarker Prediction Pipeline (CAMPP)—A standardized framework for the analysis of quantitative biological data." PLoS computational biology 16.3 (2020): e1007665.

2. **UK Biobank:**
   Sudlow, C., et al. "UK biobank: an open access resource for identifying the causes of a wide range of complex diseases of middle and old age." PLoS medicine 12.3 (2015): e1001779.

3. **Olink Proteomics:**
   Assarsson, E., et al. "Homogenous 96-plex PEA immunoassay exhibiting high sensitivity, specificity, and excellent scalability." PLoS one 9.4 (2014): e95192.

---

## Contact and Support

- **CAMPP Authors:** thilde.terkelsen@sund.ku.dk, elenap@cancer.dk
- **UK Biobank:** https://www.ukbiobank.ac.uk/
- **GitHub Issues:** https://github.com/ELELAB/CAncer-bioMarker-Prediction-Pipeline-CAMPP/issues

---

## License

This workflow and documentation are provided under the same license as CAMPP. See LICENSE.md for details.

UK Biobank data access is governed by the UK Biobank Access Management System.
