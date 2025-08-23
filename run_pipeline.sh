#!/bin/bash

# run_smn_pipeline_optimized.sh - Optimized SMN CNV Detection Pipeline
# Focus on SMN1/SMN2 exons 7&8 with MLPA-driven adaptive thresholding
# Usage: ./run_smn_pipeline_optimized.sh <bam_dir> [OPTIONS]

set -euo pipefail

# Pipeline version and metadata
PIPELINE_VERSION="2.0-SMA-OPTIMIZED"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

# Default configuration
DEFAULT_CONFIG_DIR="$PIPELINE_DIR/config"
DEFAULT_RESULTS_DIR="$PIPELINE_DIR/results_optimized"
DEFAULT_BIN_DIR="$PIPELINE_DIR/bin"
BED_FILE="$DEFAULT_CONFIG_DIR/smn_exons_critical.bed"  # Updated for exons 7&8 only
SNP_FILE="$DEFAULT_CONFIG_DIR/discriminating_snps.txt"
MLPA_TRAINING_FILE="$DEFAULT_CONFIG_DIR/mlpa_training_template.txt"

# Pipeline parameters
INPUT_BAM_DIR=""
CONFIG_DIR="$DEFAULT_CONFIG_DIR"
RESULTS_DIR="$DEFAULT_RESULTS_DIR"
BIN_DIR="$DEFAULT_BIN_DIR"
SAMPLE_TYPE="auto"
SKIP_PLOTS=false
VERBOSE=false
FORCE_THRESHOLD_RETRAIN=false
POPULATION_UPDATE=false

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m'

# Logging functions
log_header() {
    echo -e "${PURPLE}[$(date '+%Y-%m-%d %H:%M:%S')] ========================================${NC}"
    echo -e "${PURPLE}[$(date '+%Y-%m-%d %H:%M:%S')] $1${NC}"
    echo -e "${PURPLE}[$(date '+%Y-%m-%d %H:%M:%S')] ========================================${NC}"
}

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

log_step() {
    echo -e "${CYAN}[$(date '+%Y-%m-%d %H:%M:%S')] STEP: $1${NC}"
}

# Show usage information
show_usage() {
    cat << EOF
SMN CNV Detection Pipeline v${PIPELINE_VERSION}
Optimized for SMA (Spinal Muscular Atrophy) Detection

DESCRIPTION:
    This pipeline performs targeted CNV analysis of SMN1/SMN2 exons 7 and 8,
    the critical loci for SMA diagnosis and carrier screening. Features:
    
    • Adaptive MLPA-trained thresholds with continuous learning
    • Population-based clinical evidence integration
    • Fallback analysis for samples with insufficient primary coverage
    • Comprehensive SMA-specific clinical interpretation
    • Evidence-rich reporting with population frequencies

USAGE:
    $0 <bam_directory> [OPTIONS]

REQUIRED:
    bam_directory           Directory containing indexed BAM files

OPTIONS:
    --config DIR            Configuration directory (default: $DEFAULT_CONFIG_DIR)
    --results DIR           Results output directory (default: $DEFAULT_RESULTS_DIR)
    --mlpa-file FILE        MLPA training dataset (default: config/mlpa_training_template.txt)
    --sample-type TYPE      Sample type: reference|test|auto (default: auto)
    --retrain-thresholds    Force threshold retraining even if recent version exists
    --update-population     Update population evidence cache (weekly recommended)
    --skip-plots            Skip visualization generation for faster processing
    --verbose               Enable detailed logging
    --help                  Show this help message

SAMPLE TYPE AUTO-DETECTION:
    • Files containing 'ref', 'control', 'normal' → reference samples
    • All other files → test samples
    • Minimum 3 reference samples recommended for reliable normalization

KEY FEATURES:
    🧬 SMN1/SMN2 Exons 7&8 Focus: Critical loci for SMA diagnosis
    📊 MLPA-Driven Thresholds: Adaptive learning from gold-standard data
    🌍 Population Evidence: ClinVar, gnomAD, literature integration
    🔄 Fallback Analysis: Maintains interpretability with low coverage
    📋 Clinical Reports: SMA-specific interpretation with evidence context
    ⚡ Continuous Learning: Threshold optimization with each run

EXAMPLES:
    # Basic analysis with auto-detection
    $0 /data/sma_cohort/bams/
    
    # Reference cohort analysis with threshold retraining
    $0 /data/reference_samples/ --sample-type reference --retrain-thresholds
    
    # Clinical samples with population evidence update
    $0 /data/patient_samples/ --update-population --results /clinical/reports/
    
    # Fast screening mode
    $0 /data/screening_cohort/ --skip-plots --sample-type test

OUTPUT STRUCTURE:
    results_optimized/
    ├── thresholds/                 # Adaptive threshold versions with metadata
    ├── population_cache/           # Cached population evidence (ClinVar, etc.)
    ├── depth/                      # SMN exon depth extraction results
    ├── normalized/                 # Z-scores with MLPA-optimized thresholds
    ├── cnv_calls/                  # Copy number estimates with confidence
    ├── reports/                    # Per-sample TXT/HTML reports
    │   ├── SAMPLE_ID/
    │   │   ├── SAMPLE_ID_report.txt      # Detailed text report
    │   │   ├── SAMPLE_ID_report.html     # Rich HTML report with evidence
    │   │   └── SAMPLE_ID_data.json       # Structured analysis data
    ├── batch_summary_TIMESTAMP.tsv       # Consolidated batch results
    ├── batch_summary_TIMESTAMP.html      # Interactive batch dashboard
    └── logs/                       # Detailed execution logs

CLINICAL INTERPRETATION LEVELS:
    • AFFECTED: SMN1 homozygous deletion (likely SMA)
    • CARRIER: SMN1 heterozygous deletion (SMA carrier)
    • NORMAL: Standard SMN1/SMN2 copy numbers
    • UNCERTAIN: Atypical patterns requiring confirmation

QUALITY CONTROL:
    • PASS: Adequate coverage across all critical exons
    • WARNING: Some quality concerns, results interpretable
    • FAIL: Poor quality, confirmatory testing essential

For support and documentation: https://github.com/your-org/smn-cnv-pipeline
EOF
}

# Validate dependencies and configuration
validate_environment() {
    log_step "Validating environment and dependencies"
    
    local missing_tools=()
    
    # Check required tools
    for tool in samtools python3; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        exit 1
    fi
    
    # Check Python packages
    log_info "Checking Python environment..."
    python3 -c "
import sys
required_packages = ['pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'sklearn']
missing = []
for pkg in required_packages:
    try:
        __import__(pkg)
    except ImportError:
        missing.append(pkg)
if missing:
    print(f'Missing Python packages: {missing}')
    sys.exit(1)
print('All required Python packages found')
" || {
        log_error "Python package validation failed"
        echo "Install missing packages: pip install pandas numpy matplotlib seaborn scipy scikit-learn"
        exit 1
    }
    
    # Validate configuration files
    if [ ! -f "$BED_FILE" ]; then
        log_error "BED file not found: $BED_FILE"
        exit 1
    fi
    
    if [ ! -f "$SNP_FILE" ]; then
        log_error "SNP configuration file not found: $SNP_FILE"
        exit 1
    fi
    
    # Create critical exons BED file if it doesn't exist
    if [ ! -f "$BED_FILE" ]; then
        log_info "Creating optimized BED file for SMN exons 7&8..."
        mkdir -p "$(dirname "$BED_FILE")"
        cat > "$BED_FILE" << 'EOF'
# SMN1/SMN2 Critical Exons for SMA Analysis (GRCh38)
# Focus on exons 7 and 8 - essential for SMA diagnosis
chr5	70946066	70946176	SMN1_exon7	.	+
chr5	70951941	70951994	SMN1_exon8	.	+
chr5	70070641	70070751	SMN2_exon7	.	+
chr5	70076521	70076574	SMN2_exon8	.	+
EOF
    fi
    
    # Validate input BAM directory
    if [ ! -d "$INPUT_BAM_DIR" ]; then
        log_error "Input BAM directory not found: $INPUT_BAM_DIR"
        exit 1
    fi
    
    local bam_count=$(find "$INPUT_BAM_DIR" -name "*.bam" | wc -l)
    if [ "$bam_count" -eq 0 ]; then
        log_error "No BAM files found in directory: $INPUT_BAM_DIR"
        exit 1
    fi
    
    log_success "Environment validation passed ($bam_count BAM files found)"
}

# Setup output directories and initialize components
initialize_pipeline() {
    log_step "Initializing pipeline directories and components"
    
    # Create directory structure
    mkdir -p "$RESULTS_DIR"/{depth,normalized,cnv_calls,reports,logs}
    mkdir -p "$RESULTS_DIR"/{thresholds,population_cache,plots}
    
    # Initialize logging
    export PIPELINE_LOG_DIR="$RESULTS_DIR/logs"
    export PIPELINE_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
    
    log_info "Results directory: $RESULTS_DIR"
    log_info "Pipeline timestamp: $PIPELINE_TIMESTAMP"
    
    # Update population evidence if requested
    if [ "$POPULATION_UPDATE" = true ]; then
        log_info "Updating population evidence cache..."
        python3 "$BIN_DIR/update_population_evidence.py" \
            --cache-dir "$RESULTS_DIR/population_cache" \
            --force-update 2>&1 | tee "$PIPELINE_LOG_DIR/population_update.log"
    fi
    
    log_success "Pipeline initialization complete"
}

# Extract read depth for SMN critical exons
extract_smn_depth() {
    log_step "Extracting read depth for SMN critical exons 7&8"
    
    local depth_dir="$RESULTS_DIR/depth"
    local log_file="$PIPELINE_LOG_DIR/depth_extraction.log"
    
    # Use optimized depth extraction focusing on critical exons
    python3 "$BIN_DIR/smn_depth_extractor.py" \
        --bam-dir "$INPUT_BAM_DIR" \
        --bed-file "$BED_FILE" \
        --output-dir "$depth_dir" \
        --sample-type "$SAMPLE_TYPE" \
        --enable-fallback \
        --min-depth 10 \
        --min-mapq 20 \
        --min-baseq 20 \
        2>&1 | tee "$log_file"
    
    # Verify output
    local depth_files=$(find "$depth_dir" -name "*_depth_results.json" | wc -l)
    if [ "$depth_files" -eq 0 ]; then
        log_error "No depth files were created"
        exit 1
    fi
    
    log_success "Depth extraction completed ($depth_files samples processed)"
}

# Normalize coverage with MLPA-optimized thresholds
normalize_coverage_mlpa() {
    log_step "Normalizing coverage with MLPA-optimized thresholds"
    
    local log_file="$PIPELINE_LOG_DIR/normalization.log"
    
    # Load or initialize MLPA threshold manager
    local threshold_args=""
    if [ -f "$MLPA_TRAINING_FILE" ]; then
        threshold_args="--mlpa-file $MLPA_TRAINING_FILE"
    fi
    
    if [ "$FORCE_THRESHOLD_RETRAIN" = true ]; then
        threshold_args="$threshold_args --force-retrain"
    fi
    
    python3 "$BIN_DIR/mlpa_threshold_normalizer.py" \
        --depth-dir "$RESULTS_DIR/depth" \
        --output-dir "$RESULTS_DIR/normalized" \
        --threshold-dir "$RESULTS_DIR/thresholds" \
        $threshold_args \
        --population-cache "$RESULTS_DIR/population_cache" \
        2>&1 | tee "$log_file"
    
    # Verify normalization output
    if [ ! -f "$RESULTS_DIR/normalized/z_scores_optimized.txt" ]; then
        log_error "Normalization failed - Z-scores file not created"
        exit 1
    fi
    
    log_success "Coverage normalization completed with MLPA-optimized thresholds"
}

# Estimate copy numbers with adaptive thresholds
estimate_copy_numbers() {
    log_step "Estimating copy numbers with adaptive thresholds"
    
    local log_file="$PIPELINE_LOG_DIR/copy_number_estimation.log"
    
    local plot_args=""
    if [ "$SKIP_PLOTS" = false ]; then
        plot_args="--generate-plots"
    fi
    
    python3 "$BIN_DIR/adaptive_copy_number_caller.py" \
        --z-scores "$RESULTS_DIR/normalized/z_scores_optimized.txt" \
        --threshold-dir "$RESULTS_DIR/thresholds" \
        --output-dir "$RESULTS_DIR/cnv_calls" \
        --focus-exons SMN1_exon7,SMN1_exon8,SMN2_exon7,SMN2_exon8 \
        $plot_args \
        --sma-specific \
        2>&1 | tee "$log_file"
    
    # Verify copy number calling output
    if [ ! -f "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" ]; then
        log_error "Copy number estimation failed"
        exit 1
    fi
    
    log_success "Copy number estimation completed"
}

# Generate comprehensive reports with population evidence
generate_enhanced_reports() {
    log_step "Generating enhanced reports with population evidence"
    
    local log_file="$PIPELINE_LOG_DIR/report_generation.log"
    
    python3 "$BIN_DIR/enhanced_smn_reporter.py" \
        --cnv-results "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" \
        --depth-results "$RESULTS_DIR/depth" \
        --population-cache "$RESULTS_DIR/population_cache" \
        --threshold-history "$RESULTS_DIR/thresholds" \
        --output-dir "$RESULTS_DIR/reports" \
        --format all \
        --include-population-evidence \
        --sma-clinical-context \
        2>&1 | tee "$log_file"
    
    # Generate batch summary
    python3 "$BIN_DIR/batch_summary_generator.py" \
        --reports-dir "$RESULTS_DIR/reports" \
        --output-dir "$RESULTS_DIR" \
        --include-mlpa-concordance \
        --threshold-performance-analysis \
        2>&1 | tee -a "$log_file"
    
    # Count generated reports
    local report_count=$(find "$RESULTS_DIR/reports" -name "*_report.html" | wc -l)
    local batch_reports=$(find "$RESULTS_DIR" -name "batch_summary_*.html" | wc -l)
    
    log_success "Report generation completed ($report_count individual + $batch_reports batch reports)"
}

# Validate and assess pipeline results
validate_results() {
    log_step "Validating pipeline results and performing quality assessment"
    
    local validation_log="$PIPELINE_LOG_DIR/validation.log"
    
    python3 "$BIN_DIR/pipeline_validator.py" \
        --results-dir "$RESULTS_DIR" \
        --expected-samples $(find "$INPUT_BAM_DIR" -name "*.bam" | wc -l) \
        --critical-exons SMN1_exon7,SMN1_exon8,SMN2_exon7,SMN2_exon8 \
        --output-log "$validation_log" \
        2>&1 | tee "$validation_log"
    
    # Check for critical findings that need immediate attention
    local affected_samples=$(grep -c "AFFECTED" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
    local carrier_samples=$(grep -c "CARRIER" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
    local failed_samples=$(grep -c "FAIL" "$RESULTS_DIR/reports"/*/*.txt 2>/dev/null || echo "0")
    
    if [ "$affected_samples" -gt 0 ]; then
        log_warning "CRITICAL: $affected_samples samples with potential SMA-affected status detected"
    fi
    
    if [ "$carrier_samples" -gt 0 ]; then
        log_info "$carrier_samples SMA carrier samples detected"
    fi
    
    if [ "$failed_samples" -gt 0 ]; then
        log_warning "$failed_samples samples failed quality control"
    fi
    
    log_success "Results validation completed"
}

# Generate final pipeline summary
create_pipeline_summary() {
    log_step "Creating comprehensive pipeline summary"
    
    local summary_file="$RESULTS_DIR/PIPELINE_SUMMARY_${PIPELINE_TIMESTAMP}.txt"
    local start_time_file="$RESULTS_DIR/.pipeline_start_time"
    
    local runtime="Unknown"
    if [ -f "$start_time_file" ]; then
        local start_time=$(cat "$start_time_file")
        local end_time=$(date +%s)
        runtime="$((end_time - start_time)) seconds"
    fi
    
    cat > "$summary_file" << EOF
SMN CNV Detection Pipeline Summary
==================================
Pipeline Version: ${PIPELINE_VERSION}
Analysis Timestamp: $(date '+%Y-%m-%d %H:%M:%S')
Runtime: $runtime

INPUT CONFIGURATION
===================
BAM Directory: $INPUT_BAM_DIR
Configuration: $CONFIG_DIR
Results Directory: $RESULTS_DIR
Sample Type: $SAMPLE_TYPE
MLPA Training: $([ -f "$MLPA_TRAINING_FILE" ] && echo "Available" || echo "Not available")

ANALYSIS FOCUS
==============
Target Regions: SMN1/SMN2 exons 7 and 8 (critical SMA loci)
Normalization: MLPA-trained adaptive thresholds
Copy Number Calling: Population-validated with evidence integration
Clinical Context: SMA-specific interpretation with literature support

SAMPLE PROCESSING
=================
Total BAM Files: $(find "$INPUT_BAM_DIR" -name "*.bam" | wc -l)
Successfully Processed: $(find "$RESULTS_DIR/reports" -name "*_report.html" | wc -l)
Failed Quality Control: $(grep -c "FAIL" "$RESULTS_DIR/reports"/*/*.txt 2>/dev/null || echo "0")

CLINICAL FINDINGS
=================
Potential SMA Affected: $(grep -c "AFFECTED" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
SMA Carriers: $(grep -c "CARRIER" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
Normal Results: $(grep -c "NORMAL" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
Uncertain/Atypical: $(grep -c "UNCERTAIN" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")

THRESHOLD OPTIMIZATION
======================
$(if [ -f "$RESULTS_DIR/thresholds/threshold_performance.json" ]; then
    python3 -c "
import json
with open('$RESULTS_DIR/thresholds/threshold_performance.json') as f:
    perf = json.load(f)
    print(f\"Current Threshold Version: {perf.get('version', 'N/A')}\")
    print(f\"MLPA Concordance: {perf.get('mlpa_concordance', {}).get('overall', 'N/A')}\")
    print(f\"Matthews Correlation: {perf.get('mcc', 'N/A')}\")
    print(f\"Training Samples: {perf.get('training_samples', 'N/A')}\")
"
else
    echo "Threshold performance data not available"
fi)

POPULATION EVIDENCE
===================
ClinVar Integration: $([ -f "$RESULTS_DIR/population_cache/clinvar_SMN1.json" ] && echo "Updated" || echo "Not available")
Population Frequencies: $([ -f "$RESULTS_DIR/population_cache/population_frequencies.json" ] && echo "Available" || echo "Not available")
Literature Links: Integrated in individual reports

OUTPUT FILES
============
Individual Reports: $RESULTS_DIR/reports/SAMPLE_ID/
  - SAMPLE_ID_report.txt (detailed text report)
  - SAMPLE_ID_report.html (interactive HTML with evidence)
  - SAMPLE_ID_data.json (structured analysis data)

Batch Reports: $RESULTS_DIR/
  - batch_summary_${PIPELINE_TIMESTAMP}.tsv (tabular summary)
  - batch_summary_${PIPELINE_TIMESTAMP}.html (interactive dashboard)

Analysis Data:
  - $RESULTS_DIR/cnv_calls/smn_copy_numbers.txt (primary CNV calls)
  - $RESULTS_DIR/normalized/z_scores_optimized.txt (normalized data)
  - $RESULTS_DIR/thresholds/ (adaptive threshold versions)

RECOMMENDATIONS FOR NEXT STEPS
===============================
$(if [ "$(grep -c "AFFECTED" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")" -gt 0 ]; then
    echo "1. URGENT: Review samples with AFFECTED status for immediate clinical action"
fi)
$(if [ "$(grep -c "CARRIER" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")" -gt 0 ]; then
    echo "2. IMPORTANT: Genetic counseling recommended for carrier samples"
fi)
$(if [ "$(grep -c "FAIL" "$RESULTS_DIR/reports"/*/*.txt 2>/dev/null || echo "0")" -gt 0 ]; then
    echo "3. QUALITY: Consider repeat testing for failed quality control samples"
fi)
4. CONFIRMATION: MLPA or qPCR confirmation recommended for actionable findings
5. FAMILY: Consider cascade testing for at-risk family members

CLINICAL DISCLAIMERS
====================
- This analysis is for research and clinical decision support
- Confirmatory testing with validated clinical methods is recommended
- Genetic counseling should be provided for actionable findings
- Clinical correlation is essential for interpretation

For detailed analysis results, see individual sample reports in:
$RESULTS_DIR/reports/

For interactive batch overview, open:
$RESULTS_DIR/batch_summary_${PIPELINE_TIMESTAMP}.html
EOF

    log_success "Pipeline summary created: $summary_file"
    
    # Display key findings
    echo ""
    log_header "PIPELINE EXECUTION COMPLETED"
    echo ""
    log_info "Runtime: $runtime"
    log_info "Total samples processed: $(find "$RESULTS_DIR/reports" -name "*_report.html" | wc -l)"
    
    local affected=$(grep -c "AFFECTED" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
    local carriers=$(grep -c "CARRIER" "$RESULTS_DIR/cnv_calls/smn_copy_numbers.txt" 2>/dev/null || echo "0")
    
    if [ "$affected" -gt 0 ]; then
        log_warning "CRITICAL FINDINGS: $affected potential SMA-affected samples"
    fi
    
    if [ "$carriers" -gt 0 ]; then
        log_info "CLINICAL NOTE: $carriers SMA carrier samples detected"
    fi
    
    echo ""
    log_success "Results available in: $RESULTS_DIR"
    log_info "View batch summary: $RESULTS_DIR/batch_summary_${PIPELINE_TIMESTAMP}.html"
    log_info "Individual reports: $RESULTS_DIR/reports/SAMPLE_ID/"
    echo ""
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h)
                show_usage
                exit 0
                ;;
            --config)
                CONFIG_DIR="$2"
                BED_FILE="$CONFIG_DIR/smn_exons_critical.bed"
                SNP_FILE="$CONFIG_DIR/discriminating_snps.txt"
                MLPA_TRAINING_FILE="$CONFIG_DIR/mlpa_training_template.txt"
                shift 2
                ;;
            --results)
                RESULTS_DIR="$2"
                shift 2
                ;;
            --mlpa-file)
                MLPA_TRAINING_FILE="$2"
                shift 2
                ;;
            --sample-type)
                SAMPLE_TYPE="$2"
                if [[ ! "$SAMPLE_TYPE" =~ ^(reference|test|auto)$ ]]; then
                    log_error "Invalid sample type: $SAMPLE_TYPE"
                    exit 1
                fi
                shift 2
                ;;
            --retrain-thresholds)
                FORCE_THRESHOLD_RETRAIN=true
                shift
                ;;
            --update-population)
                POPULATION_UPDATE=true
                shift
                ;;
            --skip-plots)
                SKIP_PLOTS=true
                shift
                ;;
            --verbose)
                VERBOSE=true
                set -x
                shift
                ;;
            -*)
                log_error "Unknown option: $1"
                show_usage
                exit 1
                ;;
            *)
                if [ -z "$INPUT_BAM_DIR" ]; then
                    INPUT_BAM_DIR="$1"
                else
                    log_error "Multiple input directories specified"
                    exit 1
                fi
                shift
                ;;
        esac
    done
    
    if [ -z "$INPUT_BAM_DIR" ]; then
        log_error "Input BAM directory is required"
        show_usage
        exit 1
    fi
}

# Main pipeline execution
main() {
    # Record start time
    echo "$(date +%s)" > "$RESULTS_DIR/.pipeline_start_time"
    
    log_header "SMN CNV Detection Pipeline v${PIPELINE_VERSION}"
    log_info "Optimized for SMA (Spinal Muscular Atrophy) Analysis"
    log_info "Focus: SMN1/SMN2 exons 7&8 with population evidence integration"
    echo ""
    
    log_info "Input BAM directory: $INPUT_BAM_DIR"
    log_info "Results directory: $RESULTS_DIR"
    log_info "Sample type: $SAMPLE_TYPE"
    log_info "MLPA training file: $([ -f "$MLPA_TRAINING_FILE" ] && echo "$MLPA_TRAINING_FILE" || echo "Not provided")"
    echo ""
    
    # Execute pipeline steps
    validate_environment
    initialize_pipeline
    extract_smn_depth
    normalize_coverage_mlpa  
    estimate_copy_numbers
    generate_enhanced_reports
    validate_results
    create_pipeline_summary
    
    # Cleanup
    rm -f "$RESULTS_DIR/.pipeline_start_time"
    
    log_success "SMN CNV Detection Pipeline completed successfully!"
}

# Script entry point
parse_arguments "$@"
main
