# Create discriminating SNPs file
    cat > "$PIPELINE_ROOT/config/discriminating_snps.txt" << 'EOF'
# Discriminating SNPs for SMN1/SMN2 differentiation
# Format: chrom	pos	ref	alt	gene	exon	type
chr5	70946146	C	T	SMN1	exon7	discriminating
chr5	70946169	G	A	SMN1	exon7	discriminating
chr5	70951986	C	T	SMN1	exon8	discriminating
chr5	70070721	C	T	SMN2	exon7	discriminating
chr5	70070744	G	A	SMN2	exon7	discriminating
chr5	70076566	C	T	SMN2	exon8	discriminating
EOF

    # Create MLPA training template
    cat > "$PIPELINE_ROOT/config/mlpa_training_template.txt" << 'EOF'
# MLPA Training Dataset Template for SMN CNV Analysis
# Format: sample_id	SMN1_copy_number	SMN2_copy_number	clinical_interpretation	validation_method
CTRL001	2	2	NORMAL	MLPA_validated
CTRL002	2	3	NORMAL	MLPA_validated
CTRL003	2	2	NORMAL	MLPA_validated
CTRL004	2	1	NORMAL	MLPA_validated
CTRL005	2	4	NORMAL	MLPA_validated
CARRIER001	1	2	CARRIER	MLPA_validated
CARRIER002	1	3	CARRIER	MLPA_validated
CARRIER003	1	2	CARRIER	MLPA_validated
CARRIER004	1	1	CARRIER	MLPA_validated
CARRIER005	1	4	CARRIER	MLPA_validated
AFFECTED001	0	2	AFFECTED	MLPA_validated
AFFECTED002	0	3	AFFECTED	MLPA_validated
AFFECTED003	0	1	AFFECTED_TYPE1	MLPA_validated
AFFECTED004	0	4	AFFECTED_TYPE4	MLPA_validated
AFFECTED005	0	2	AFFECTED	MLPA_validated
EOF

    # Create pipeline configuration
    cat > "$PIPELINE_ROOT/config/pipeline_config.json" << 'EOF'
{
    "pipeline_version": "2.0-SMA-OPTIMIZED",
    "target_genes": ["SMN1", "SMN2"],
    "critical_exons": ["SMN1_exon7", "SMN1_exon8", "SMN2_exon7", "SMN2_exon8"],
    "reference_genome": "GRCh38",
    "quality_thresholds": {
        "min_depth": 10,
        "min_mapq": 20,
        "min_baseq": 20,
        "min_coverage_fraction": 0.8,
        "min_quality_score": 0.5
    },
    "cnv_thresholds": {
        "deletion_z_threshold": -1.5,
        "duplication_z_threshold": 1.5,
        "confidence_threshold": 0.7
    },
    "sma_interpretation": {
        "enable_sma_specific": true,
        "severity_prediction": true,
        "population_evidence": true
    }
}
EOF
    
    log_success "Configuration files created"
}

# Make Python scripts executable and set up bin directory
setup_bin_directory() {
    log_info "Setting up bin directory with Python scripts"
    
    # Note: In real setup, you would copy the actual Python files here
    # For this example, we'll create placeholder files with instructions
    
    cat > "$PIPELINE_ROOT/bin/README.md" << 'EOF'
# SMN Pipeline Python Scripts

This directory should contain the following Python scripts:

## Required Scripts:
1. `smn_depth_extractor.py` - Extract read depth for SMN exons
2. `mlpa_threshold_normalizer.py` - MLPA-optimized threshold normalization  
3. `adaptive_copy_number_caller.py` - Adaptive copy number calling
4. `enhanced_smn_reporter.py` - Generate comprehensive reports
5. `batch_summary_generator.py` - Create batch summaries
6. `pipeline_validator.py` - Validate pipeline results
7. `update_population_evidence.py` - Update population databases

## Setup Instructions:
1. Copy the Python scripts provided by the assistant to this directory
2. Make them executable: `chmod +x *.py`
3. Ensure all dependencies are installed (see requirements.txt)
4. Test with: `python3 script_name.py --help`

## Directory Structure:
```
bin/
├── smn_depth_extractor.py
├── mlpa_threshold_normalizer.py
├── adaptive_copy_number_caller.py
├── enhanced_smn_reporter.py
├── batch_summary_generator.py
├── pipeline_validator.py
├── update_population_evidence.py
└── README.md
```
EOF

    log_warning "Python scripts need to be manually copied to $PIPELINE_ROOT/bin/"
    log_info "See $PIPELINE_ROOT/bin/README.md for instructions"
}

# Create test data and documentation
setup_documentation() {
    log_info "Setting up documentation and test data"
    
    # Create main README
    cat > "$PIPELINE_ROOT/README.md" << 'EOF'
# SMN CNV Detection Pipeline v2.0

## Overview
This pipeline performs targeted CNV analysis of SMN1/SMN2 exons 7 and 8 for SMA (Spinal Muscular Atrophy) detection and carrier screening.

## Key Features
- 🧬 **SMN1/SMN2 Focus**: Critical exons 7&8 analysis
- 📊 **MLPA Integration**: Adaptive thresholds from gold-standard data
- 🌍 **Population Evidence**: ClinVar and literature integration
- 🔄 **Adaptive Learning**: Continuous threshold optimization
- 📋 **Clinical Reports**: SMA-specific interpretation

## Quick Start
```bash
# Setup pipeline
./setup_pipeline.sh

# Run analysis
./run_smn_pipeline_optimized.sh /path/to/bam/files/

# View results
open results_optimized/batch_summary_TIMESTAMP.html
```

## Requirements
- samtools
- Python 3.7+
- Required Python packages (see requirements.txt)
- Indexed BAM files
- MLPA training data (optional but recommended)

## Directory Structure
```
cnv_pipeline_ml/
├── bin/                    # Python scripts
├── config/                 # Configuration files
├── results_optimized/      # Analysis results
├── docs/                   # Documentation
├── tests/                  # Test data and scripts
├── run_smn_pipeline_optimized.sh  # Main pipeline script
└── setup_pipeline.sh       # Setup script
```

## Output
- Individual sample reports (TXT/HTML/JSON)
- Batch summary dashboard (HTML)
- Copy number calls (TSV)
- Quality metrics and validation logs

## Clinical Use
This pipeline is designed for research and clinical decision support. All actionable findings should be confirmed with validated clinical testing methods (MLPA, qPCR).

## Support
For issues and questions, please refer to the documentation or contact the development team.
EOF

    # Create example test script
    cat > "$PIPELINE_ROOT/tests/test_pipeline.sh" << 'EOF'
#!/bin/bash
# Test script for SMN pipeline
# Creates minimal test data and runs pipeline validation

echo "SMN Pipeline Test Suite"
echo "======================"

# Create test BAM directory with mock files (you'll need real BAM files)
mkdir -p test_data/bams
echo "Create test BAM files in test_data/bams/ for testing"

# Test configuration validation
echo "Testing configuration files..."
if [ -f "../config/smn_exons.bed" ]; then
    echo "✓ BED file found"
else
    echo "✗ BED file missing"
fi

if [ -f "../config/discriminating_snps.txt" ]; then
    echo "✓ SNP file found"
else
    echo "✗ SNP file missing"
fi

echo "Test script complete. Add real BAM files to test_data/bams/ to run full pipeline test."
EOF

    chmod +x "$PIPELINE_ROOT/tests/test_pipeline.sh"
    
    log_success "Documentation and test framework created"
}

# Validate installation
validate_installation() {
    log_info "Validating pipeline installation"
    
    local issues=()
    
    # Check directory structure
    for dir in bin config results_optimized docs tests; do
        if [ ! -d "$PIPELINE_ROOT/$dir" ]; then
            issues+=("Missing directory: $dir")
        fi
    done
    
    # Check main script
    if [ ! -f "$PIPELINE_ROOT/run_smn_pipeline_optimized.sh" ]; then
        issues+=("Missing main pipeline script")
    fi
    
    # Check configuration files
    for file in config/smn_exons.bed config/discriminating_snps.txt config/mlpa_training_template.txt; do
        if [ ! -f "$PIPELINE_ROOT/$file" ]; then
            issues+=("Missing config file: $file")
        fi
    done
    
    # Check Python environment
    if ! python3 -c "import pandas, numpy, matplotlib, seaborn, scipy, sklearn" 2>/dev/null; then
        issues+=("Python dependencies not properly installed")
    fi
    
    if [ ${#issues[@]} -eq 0 ]; then
        log_success "Pipeline installation validated successfully"
        return 0
    else
        log_error "Installation validation failed:"
        for issue in "${issues[@]}"; do
            echo "  - $issue"
        done
        return 1
    fi
}

# Main setup function
main() {
    echo "=================================================="
    echo "SMN CNV Detection Pipeline Setup v2.0"
    echo "Optimized for SMA Analysis"
    echo "=================================================="
    echo ""
    
    log_info "Pipeline root directory: $PIPELINE_ROOT"
    echo ""
    
    # Create directories
    create_directories
    
    # Install dependencies
    install_system_deps
    install_python_deps
    
    # Setup configuration
    setup_config_files
    setup_bin_directory
    setup_documentation
    
    echo ""
    log_info "Setup completed. Next steps:"
    echo ""
    echo "1. Copy Python scripts to $PIPELINE_ROOT/bin/"
    echo "   - Use the Python files provided by the assistant"
    echo "   - Make them executable: chmod +x $PIPELINE_ROOT/bin/*.py"
    echo ""
    echo "2. Copy the main pipeline script:"
    echo "   - Save run_smn_pipeline_optimized.sh to $PIPELINE_ROOT/"
    echo "   - Make executable: chmod +x $PIPELINE_ROOT/run_smn_pipeline_optimized.sh"
    echo ""
    echo "3. Update MLPA training data:"
    echo "   - Edit $PIPELINE_ROOT/config/mlpa_training_template.txt"
    echo "   - Add your actual MLPA validation samples"
    echo ""
    echo "4. Test the installation:"
    echo "   - cd $PIPELINE_ROOT && ./tests/test_pipeline.sh"
    echo ""
    echo "5. Run pipeline on your data:"
    echo "   - ./run_smn_pipeline_optimized.sh /path/to/bam/files/"
    echo ""
    
    # Validate installation
    if validate_installation; then
        log_success "Setup completed successfully!"
        echo ""
        echo "Pipeline is ready for use. See README.md for detailed instructions."
    else
        log_warning "Setup completed with issues. Please address the problems above."
    fi
}

# Check if running as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi#!/bin/bash

# setup_pipeline.sh - Setup SMN CNV Detection Pipeline
# Creates directory structure and installs dependencies

set -euo pipefail

PIPELINE_ROOT="/data/Sethu/SMN/cnv_pipeline_ml"
BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1${NC}"
}

log_success() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $1${NC}"
}

log_warning() {
    echo -e "${YELLOW}[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: $1${NC}"
}

log_error() {
    echo -e "${RED}[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1${NC}"
}

# Create directory structure
create_directories() {
    log_info "Creating pipeline directory structure"
    
    mkdir -p "$PIPELINE_ROOT"/{bin,config,results_optimized,docs,tests}
    mkdir -p "$PIPELINE_ROOT/results_optimized"/{depth,normalized,cnv_calls,reports,logs}
    mkdir -p "$PIPELINE_ROOT/results_optimized"/{thresholds,population_cache,plots}
    
    log_success "Directory structure created"
}

# Install system dependencies
install_system_deps() {
    log_info "Checking system dependencies"
    
    # Check if samtools is installed
    if ! command -v samtools &> /dev/null; then
        log_warning "samtools not found. Please install samtools:"
        echo "  Ubuntu/Debian: sudo apt-get install samtools"
        echo "  CentOS/RHEL: sudo yum install samtools"
        echo "  macOS: brew install samtools"
        echo "  Conda: conda install -c bioconda samtools"
    else
        log_success "samtools found: $(samtools --version | head -n1)"
    fi
    
    # Check Python version
    if ! command -v python3 &> /dev/null; then
        log_error "Python 3 is required but not found"
        exit 1
    else
        python_version=$(python3 --version)
        log_success "Python found: $python_version"
    fi
}

# Install Python dependencies
install_python_deps() {
    log_info "Installing Python dependencies"
    
    # Create requirements.txt
    cat > "$PIPELINE_ROOT/requirements.txt" << 'EOF'
pandas>=1.3.0
numpy>=1.21.0
matplotlib>=3.5.0
seaborn>=0.11.0
scipy>=1.7.0
scikit-learn>=1.0.0
plotly>=5.0.0
jinja2>=3.0.0
requests>=2.25.0
biopython>=1.79
pysam>=0.19.0
EOF

    # Install packages
    log_info "Installing Python packages..."
    python3 -m pip install -r "$PIPELINE_ROOT/requirements.txt"
    
    log_success "Python dependencies installed"
}

# Create configuration files
setup_config_files() {
    log_info "Setting up configuration files"
    
    # Create smn_exons.bed
    cat > "$PIPELINE_ROOT/config/smn_exons.bed" << 'EOF'
# SMN1/SMN2 Critical Exons for SMA Analysis (GRCh38)
# Focus on exons 7 and 8 - essential for SMA diagnosis
chr5	70946066	70946176	SMN1_exon7	.	+
chr5	70951941	70951994	SMN1_exon8	.	+
chr5	70070641	70070751	SMN2_exon7	.	+
chr5	70076521	70076574	SMN2_exon8	.	+
EOF

    # Create discriminating SNPs file
    cat > "$PIPELINE_ROOT/config/discriminating_snps.txt" << 'EOF'
# Discriminating SNPs for SMN1/SMN2 differentiation
# Format: chrom	pos	ref	alt	gene	exon	type
chr5	70
