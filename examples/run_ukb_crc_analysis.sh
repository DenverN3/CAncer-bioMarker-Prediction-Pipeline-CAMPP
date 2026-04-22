#!/bin/bash
#
# UK Biobank Colorectal Cancer Biomarker Discovery - Complete Analysis Pipeline
#
# This script runs the complete CAMPP workflow for UK Biobank CRC proteomics data
# Author: CAMPP Team
# Date: 2025-12-06
#
# Prerequisites:
# 1. R version 4.0.0 or higher installed
# 2. Input files: ukb_crc_proteomics.txt and ukb_crc_metadata.txt
# 3. CAMPP.R and CAMPPFunctions.R in the same directory
#
# Usage:
#   bash run_ukb_crc_analysis.sh
#
# All results will be saved in the Results/ directory

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Project name
PROJECT="UKB_CRC_Complete"

# Data files
DATA_FILE="ukb_crc_proteomics.txt"
META_FILE="ukb_crc_metadata.txt"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}UK Biobank CRC Analysis Pipeline${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Check if input files exist
if [ ! -f "$DATA_FILE" ]; then
    echo -e "${RED}ERROR: Data file not found: $DATA_FILE${NC}"
    echo "Please ensure the proteomics data file is in the current directory"
    echo "You can use the template: examples/ukb_crc_proteomics_template.txt"
    exit 1
fi

if [ ! -f "$META_FILE" ]; then
    echo -e "${RED}ERROR: Metadata file not found: $META_FILE${NC}"
    echo "Please ensure the metadata file is in the current directory"
    echo "You can use the template: examples/ukb_crc_metadata_template.txt"
    exit 1
fi

# Check if CAMPP.R exists
if [ ! -f "CAMPP.R" ]; then
    echo -e "${RED}ERROR: CAMPP.R not found in current directory${NC}"
    exit 1
fi

echo -e "${GREEN}✓${NC} Input files found"
echo ""

# Create results directory
mkdir -p Results
echo -e "${GREEN}✓${NC} Results directory created"
echo ""

# ============================================================================
# STEP 1: Quality Control and Data Exploration
# ============================================================================
echo -e "${YELLOW}[Step 1/6]${NC} Running Quality Control and Data Exploration..."
echo "This step will:"
echo "  - Check data distributions (Cullen-Frey plots)"
echo "  - Generate MDS plot for sample clustering"
echo "  - Apply log2 transformation and median centering"
echo "  - Correct for batch effects (Olink plates)"
echo ""

Rscript CAMPP.R \
  -d $DATA_FILE \
  -m $META_FILE \
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

echo -e "${GREEN}✓${NC} Step 1 complete"
echo ""

# ============================================================================
# STEP 2: Differential Expression Analysis
# ============================================================================
echo -e "${YELLOW}[Step 2/6]${NC} Running Differential Expression Analysis..."
echo "This step will:"
echo "  - Identify proteins differentially expressed between CRC and controls"
echo "  - Adjust for age, sex, BMI, and smoking status"
echo "  - Apply significance cutoffs: |logFC| > 0.5, FDR < 0.05"
echo "  - Generate heatmap of significant proteins"
echo ""

Rscript CAMPP.R \
  -d $DATA_FILE \
  -m $META_FILE \
  -v other \
  -g ids,group \
  -n ${PROJECT}_DE \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -r TRUE,age,sex,bmi,smoking \
  -a DE

echo -e "${GREEN}✓${NC} Step 2 complete"
echo "  → Check Results/DEAAResults/${PROJECT}_DE_DE.txt for DE results"
echo ""

# ============================================================================
# STEP 3: LASSO Regularized Regression
# ============================================================================
echo -e "${YELLOW}[Step 3/6]${NC} Running LASSO Feature Selection..."
echo "This step will:"
echo "  - Select minimal protein signature using LASSO regression"
echo "  - Perform 10-fold cross-validation"
echo "  - Calculate classification AUC"
echo "  - Identify consensus biomarkers (DE ∩ LASSO)"
echo ""

Rscript CAMPP.R \
  -d $DATA_FILE \
  -m $META_FILE \
  -v other \
  -g ids,group \
  -n ${PROJECT}_LASSO \
  -z median \
  -t log2 \
  -b batch \
  -l 1.0 \
  -a Consensus

echo -e "${GREEN}✓${NC} Step 3 complete"
echo "  → Check Results/LASSOResults/${PROJECT}_LASSO_LASSO.txt for selected proteins"
echo "  → Check Results/LASSOResults/${PROJECT}_LASSO_DEA_LASSO_Consensus.txt for consensus biomarkers"
echo ""

# ============================================================================
# STEP 4: Survival Analysis
# ============================================================================
echo -e "${YELLOW}[Step 4/6]${NC} Running Survival Analysis..."
echo "This step will:"
echo "  - Test consensus biomarkers for survival prediction"
echo "  - Calculate hazard ratios with 95% CI"
echo "  - Adjust for clinical covariates (age, sex, BMI, smoking, stage)"
echo "  - Generate forest plots"
echo ""

Rscript CAMPP.R \
  -d $DATA_FILE \
  -m $META_FILE \
  -v other \
  -g ids,group \
  -n ${PROJECT}_Survival \
  -z median \
  -t log2 \
  -b batch \
  -u Consensus \
  -r TRUE,age,sex,bmi,smoking,stage \
  -q 20

echo -e "${GREEN}✓${NC} Step 4 complete"
echo "  → Check Results/SurvivalResults/${PROJECT}_Survival_survival.txt for hazard ratios"
echo ""

# ============================================================================
# STEP 5: WGCNA Co-expression Network Analysis
# ============================================================================
echo -e "${YELLOW}[Step 5/6]${NC} Running WGCNA Network Analysis..."
echo "This step will:"
echo "  - Identify co-expressed protein modules"
echo "  - Calculate intramodular connectivity scores"
echo "  - Generate module dendrograms and heatmaps"
echo "  - Return top 30% hub proteins per module"
echo ""

Rscript CAMPP.R \
  -d $DATA_FILE \
  -m $META_FILE \
  -v other \
  -g ids,group \
  -n ${PROJECT}_WGCNA \
  -z median \
  -t log2 \
  -b batch \
  -w DE \
  -x 5,25,30

echo -e "${GREEN}✓${NC} Step 5 complete"
echo "  → Check Results/WGCNAResults/WGCNAPlots/ for network visualizations"
echo ""

# ============================================================================
# STEP 6: Protein-Protein Interaction Networks
# ============================================================================
echo -e "${YELLOW}[Step 6/6]${NC} Running PPI Network Analysis..."
echo "This step will:"
echo "  - Map DE proteins to STRING database"
echo "  - Extract protein-protein interactions"
echo "  - Generate network visualization"
echo ""

Rscript CAMPP.R \
  -d $DATA_FILE \
  -m $META_FILE \
  -v other \
  -g ids,group \
  -n ${PROJECT}_PPI \
  -z median \
  -t log2 \
  -b batch \
  -f 0.5,0.05 \
  -p hgnc_symbol,11.0

echo -e "${GREEN}✓${NC} Step 6 complete"
echo "  → Check Results/InteractionResults/ for PPI networks"
echo ""

# ============================================================================
# Summary
# ============================================================================
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Analysis Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Results have been saved to the following directories:"
echo ""
echo "  📊 Data Quality:        DataChecks/"
echo "  📈 MDS Plots:           ${PROJECT}_*_MDSplot.pdf"
echo "  🧬 Differential Expr:   Results/DEAAResults/"
echo "  🎯 LASSO Selection:     Results/LASSOResults/"
echo "  ⏱  Survival Analysis:   Results/SurvivalResults/"
echo "  🕸  WGCNA Modules:       Results/WGCNAResults/"
echo "  🔗 PPI Networks:        Results/InteractionResults/"
echo "  🔥 Heatmaps:            ${PROJECT}_*_heatmap.pdf"
echo "  📝 Execution Log:       CAMPPlog.txt"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Review CAMPPlog.txt for any warnings or errors"
echo "  2. Check MDS plots to ensure no outliers or batch effects remain"
echo "  3. Examine consensus biomarkers in:"
echo "     Results/LASSOResults/${PROJECT}_LASSO_DEA_LASSO_Consensus.txt"
echo "  4. Review survival results for prognostic biomarkers"
echo "  5. Validate findings in an independent cohort"
echo ""
echo -e "${GREEN}For detailed interpretation, see UK_BIOBANK_USE_CASE.md${NC}"
echo ""
