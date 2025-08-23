# SMN CNV Detection Pipeline Setup Guide

## Overview
This guide will help you set up and run the SMN CNV Detection Pipeline v2.0, which is optimized for Spinal Muscular Atrophy (SMA) analysis focusing on SMN1/SMN2 exons 7 and 8.

## Prerequisites

### System Requirements
- Linux/Unix environment (tested on Ubuntu/CentOS)
- Bash shell
- At least 8GB RAM (16GB recommended for large cohorts)
- 50GB+ free disk space for results

### Required Software
1. **samtools** (for BAM file processing)
   ```bash
   # Ubuntu/Debian
   sudo apt-get install samtools
   
   # CentOS/RHEL
   sudo yum install samtools
   
   # macOS
   brew install samtools
   
   # Conda
   conda install -c bioconda samtools
   ```

2. **Python 3.7+** with pip
   ```bash
   python3 --version  # Should be 3.7 or higher
   ```

## Installation Steps

### Step 1: Create Pipeline Directory Structure
```bash
# Set your desired pipeline location
export PIPELINE_ROOT="/data/Sethu/SMN/cnv_pipeline_ml"

# Create directory structure
mkdir -p $PIPELINE_ROOT/{bin,config,results_optimized,docs,tests}
mkdir -p $PIPELINE_ROOT/results_optimized/{depth,normalized,cnv_calls,reports,logs}
mkdir -p $PIPELINE_ROOT/results_optimized/{thresholds,population_cache,plots}

cd $PIPELINE_ROOT
```

### Step 2: Install Python Dependencies
```bash
# Create and activate virtual environment (recommended)
python3 -m venv smn_pipeline_env
source smn_pipeline_env/bin/activate

# Install required packages
pip install pandas numpy matplotlib seaborn scipy scikit-learn plotly jinja2 requests biopython pysam
```

### Step 3: Download and Setup Scripts

#### A. Save the main pipeline script
Save the provided `run_smn_pipeline_optimized.sh` to `$PIPELINE_ROOT/run_smn_pipeline_optimized.sh`:
```bash
chmod +x run_smn_pipeline_optimized.sh
```

#### B. Save Python scripts to bin/ directory
Copy these Python scripts to `$PIPELINE_ROOT/bin/`:
- `smn_depth_extractor.py`
- `mlpa_threshold_normalizer.py`
- `adaptive_copy_number_caller.py`
- `enhanced_smn_reporter.py`
- `batch_summary_generator.py`
- `pipeline_validator.py`
- `update_population_evidence.py`

```bash
# Make all Python scripts executable
chmod +x bin/*.py
```

#### C. Setup configuration files
Create these files in `$PIPELINE_ROOT/config/`:

**smn_exons.bed:**
```
chr5	70946066	70946176	SMN1_exon7	.	+
chr5	70951941	70951994	SMN1_exon8	.	+
chr5	70070641	70070751	SMN2_exon7	.	+
chr5	70076521	70076574	SMN2_exon8	.	+
```

**discriminating_snps.txt:**
```
chr5	70946146	C	T	SMN1	exon7	discriminating
chr5	70946169	G	A	SMN1	exon7	discriminating
chr5	70951986	C	T	SMN1	exon8	discriminating
chr5	70070721	C	T	SMN2	exon7	discriminating
chr5	70070744	G	A	SMN2	exon7	discriminating
chr5	70076566	C	T	SMN2	exon8	discriminating
```

**mlpa_training_template.txt:**
```
sample_id	SMN1_copy_number	SMN2_copy_number	clinical_interpretation	validation_method
CTRL001	2	2	NORMAL	MLPA_validated
CARRIER001	1	2	CARRIER	MLPA_validated
AFFECTED001	0	2	AFFECTED	MLPA_validated
# Add your actual MLPA data here
```

### Step 4: Prepare Your Data

#### BAM File Requirements
- **Indexed BAM files** (.bam + .bai files)
- **Mapped to GRCh38** reference genome
- **Good coverage** over SMN regions (chromosome 5)
- **Proper naming**: Use descriptive sample IDs

#### BAM File Organization
```bash
# Organize your BAM files
mkdir -p /data/your_project/bams/
# Copy your BAM files and indices here

# Example structure:
/data/your_project/bams/
├── sample001.bam
├── sample001.bam.bai
├── sample002.bam
├── sample002.bam.bai
└── ...
```

#### MLPA Training Data (Optional but Recommended)
If you have MLPA validation data, update `config/mlpa_training_template.txt`:
- Include at least 10 samples per category (NORMAL, CARRIER, AFFECTED)
- Ensure high-quality MLPA results
- Include diverse SMN2 copy numbers

## Running the Pipeline

### Basic Usage
```bash
cd $PIPELINE_ROOT

# Basic analysis with auto-detection
./run_smn_pipeline_optimized.sh /data/your_project/bams/

# With custom output directory
./run_smn_pipeline_optimized.sh /data/your_project/bams/ --results /data/results/

# Reference samples with threshold training
./run_smn_pipeline_optimized.sh /data/reference_samples/ --sample-type reference --retrain-thresholds

# Fast screening mode
./run_smn_pipeline_optimized.sh /data/screening/ --skip-plots --sample-type test
```

### Advanced Options
```bash
# Full featured run with population update
./run_smn_pipeline_optimized.sh /data/bams/ \
    --update-population \
    --retrain-thresholds \
    --results /data/sma_analysis/ \
    --verbose

# Custom MLPA file
./run_smn_pipeline_optimized.sh /data/bams/ \
    --mlpa-file /data/your_mlpa_data.txt \
    --retrain-thresholds
```

## Understanding Results

### Directory Structure After Analysis
```
results_optimized/
├── depth/                          # Raw depth extraction
├── normalized/                     # Z-score normalized data
├── cnv_calls/                      # Copy number calls
├── reports/                        # Individual sample reports
│   └── SAMPLE_ID/
│       ├── SAMPLE_ID_report.txt    # Text report
│       ├── SAMPLE_ID_report.html   # Interactive report
│       └── SAMPLE_ID_data.json     # Structured data
├── batch_summary_TIMESTAMP.html    # Batch dashboard
├── batch_summary_TIMESTAMP.tsv     # Batch summary table
└── logs/                           # Execution logs
```

### Key Output Files
1. **Batch Dashboard**: `batch_summary_TIMESTAMP.html` - Interactive overview
2. **Copy Number Calls**: `cnv_calls/smn_copy_numbers.txt` - Main results table
3. **Individual Reports**: `reports/SAMPLE_ID/` - Detailed per-sample analysis
4. **Pipeline Summary**: `PIPELINE_SUMMARY_TIMESTAMP.txt` - Execution summary

### Clinical Interpretation
- **AFFECTED**: SMN1 homozygous deletion (likely SMA)
- **CARRIER**: SMN1 heterozygous deletion (SMA carrier)  
- **NORMAL**: Standard SMN1/SMN2 copy numbers
- **UNCERTAIN**: Atypical patterns requiring confirmation

### Quality Flags
- **PASS**: Adequate coverage, reliable results
- **WARNING**: Some concerns, interpret with caution
- **FAIL**: Poor quality, confirmatory testing essential

## Troubleshooting

### Common Issues

#### 1. "No BAM files found"
- Check BAM directory path
- Ensure files have .bam extension
- Verify read permissions

#### 2. "samtools not found"
- Install samtools (see prerequisites)
- Check PATH environment variable

#### 3. "Python package missing"
- Activate virtual environment
- Install missing packages: `pip install package_name`

#### 4. "Permission denied"
- Make scripts executable: `chmod +x *.sh bin/*.py`
- Check directory permissions

#### 5. Low quality results
- Check BAM file quality and coverage
- Verify reference genome alignment (should be GRCh38)
- Review SMN region coverage specifically

### Performance Optimization
- **Large cohorts (>100 samples)**: Use `--skip-plots` for faster processing
- **Limited resources**: Process in smaller batches
- **Reference samples**: Run with `--sample-type reference` first to establish baseline

## Validation and Quality Control

### Before Clinical Use
1. **Test with known samples** (MLPA-validated controls)
2. **Validate thresholds** with your population
3. **Review concordance** with existing clinical methods
4. **Establish SOPs** for result interpretation

### Ongoing Quality Monitoring
- **Weekly**: Update population evidence (`--update-population`)
- **Monthly**: Retrain thresholds with new MLPA data
- **Quarterly**: Review threshold performance and adjust if needed

## Clinical Considerations

### Important Notes
- This pipeline is for **research and clinical decision support**
- **Confirmatory testing** required for actionable findings
- **Genetic counseling** recommended for carriers and affected individuals
- **Family screening** indicated for positive results

### Follow-up Actions
- **AFFECTED samples**: Urgent genetic counseling, confirmatory testing
- **CARRIER samples**: Genetic counseling, partner testing if reproductive planning
- **FAIL quality**: Repeat testing with fresh sample

## Support and Documentation

### Getting Help
1. Check the logs in `results_optimized/logs/`
2. Review the pipeline summary file
3. Validate installation with test scripts
4. Consult documentation in `docs/` directory

### File Locations
- **Main script**: `run_smn_pipeline_optimized.sh`
- **Python modules**: `bin/`
- **Configuration**: `config/`
- **Results**: `results_optimized/`
- **Documentation**: `docs/`
- **Tests**: `tests/`

This pipeline represents a research tool optimized for SMA analysis. Always follow institutional guidelines and confirm critical findings with validated clinical methods.
