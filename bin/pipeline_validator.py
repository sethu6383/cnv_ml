#!/usr/bin/env python3
"""
Pipeline Validator - Validate SMN pipeline results and assess quality
"""

import os
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
import logging
from datetime import datetime

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)

class PipelineValidator:
    """Validate SMN pipeline execution and results"""
    
    def __init__(self, results_dir):
        self.results_dir = Path(results_dir)
        self.validation_results = {}
        
    def validate_directory_structure(self):
        """Validate expected directory structure exists"""
        required_dirs = [
            'depth', 'normalized', 'cnv_calls', 'reports', 
            'logs', 'thresholds', 'population_cache'
        ]
        
        missing_dirs = []
        
        for dir_name in required_dirs:
            dir_path = self.results_dir / dir_name
            if not dir_path.exists():
                missing_dirs.append(dir_name)
        
        self.validation_results['directory_structure'] = {
            'status': 'PASS' if not missing_dirs else 'FAIL',
            'missing_directories': missing_dirs,
            'required_directories': required_dirs
        }
        
        return len(missing_dirs) == 0
    
    def validate_file_completeness(self, expected_samples):
        """Validate that all expected output files were created"""
        validation = {
            'depth_files': self._count_files('depth', '*_depth_results.json'),
            'reports': self._count_files('reports', '*_report.html'),
            'cnv_calls': (self.results_dir / 'cnv_calls' / 'smn_copy_numbers.txt').exists(),
            'normalized_data': (self.results_dir / 'normalized' / 'z_scores_optimized.txt').exists(),
            'thresholds': (self.results_dir / 'thresholds' / 'current_thresholds.json').exists()
        }
        
        # Check completeness
        depth_complete = validation['depth_files'] >= expected_samples * 0.9  # Allow 10% failure
        reports_complete = validation['reports'] >= expected_samples * 0.9
        
        self.validation_results['file_completeness'] = {
            'status': 'PASS' if (depth_complete and reports_complete and 
                                validation['cnv_calls'] and validation['normalized_data']) else 'FAIL',
            'expected_samples': expected_samples,
            'depth_files_created': validation['depth_files'],
            'reports_created': validation['reports'],
            'cnv_calls_file': validation['cnv_calls'],
            'normalized_data_file': validation['normalized_data'],
            'thresholds_file': validation['thresholds']
        }
        
        return self.validation_results['file_completeness']['status'] == 'PASS'
    
    def _count_files(self, subdir, pattern):
        """Count files matching pattern in subdirectory"""
        dir_path = self.results_dir / subdir
        if not dir_path.exists():
            return 0
        return len(list(dir_path.glob(pattern)))
    
    def validate_critical_exons_coverage(self, critical_exons):
        """Validate that critical exons were properly analyzed"""
        cnv_file = self.results_dir / 'cnv_calls' / 'smn_copy_numbers.txt'
        
        if not cnv_file.exists():
            self.validation_results['critical_exons'] = {
                'status': 'FAIL',
                'reason': 'CNV calls file not found'
            }
            return False
        
        try:
            cnv_df = pd.read_csv(cnv_file, sep='\t')
            
            # Check for required exon columns
            exon_list = critical_exons.split(',') if isinstance(critical_exons, str) else critical_exons
            missing_exons = []
            
            for exon in exon_list:
                z_score_col = f"{exon}_z_score"
                cn_col = f"{exon}_copy_number"
                if z_score_col not in cnv_df.columns or cn_col not in cnv_df.columns:
                    missing_exons.append(exon)
            
            # Check for valid data
            samples_with_data = 0
            for _, row in cnv_df.iterrows():
                has_valid_data = True
                for exon in exon_list:
                    z_score_col = f"{exon}_z_score"
                    if z_score_col in cnv_df.columns:
                        if pd.isna(row[z_score_col]):
                            has_valid_data = False
                            break
                if has_valid_data:
                    samples_with_data += 1
            
            data_completeness = samples_with_data / len(cnv_df) if len(cnv_df) > 0 else 0
            
            self.validation_results['critical_exons'] = {
                'status': 'PASS' if not missing_exons and data_completeness >= 0.8 else 'FAIL',
                'missing_exons': missing_exons,
                'data_completeness': round(data_completeness, 3),
                'samples_with_valid_data': samples_with_data,
                'total_samples': len(cnv_df)
            }
            
            return self.validation_results['critical_exons']['status'] == 'PASS'
            
        except Exception as e:
            self.validation_results['critical_exons'] = {
                'status': 'FAIL',
                'reason': f'Error reading CNV file: {e}'
            }
            return False
    
    def validate_quality_metrics(self):
        """Validate overall pipeline quality metrics"""
        cnv_file = self.results_dir / 'cnv_calls' / 'smn_copy_numbers.txt'
        
        if not cnv_file.exists():
            self.validation_results['quality_metrics'] = {
                'status': 'FAIL',
                'reason': 'CNV calls file not found'
            }
            return False
        
        try:
            cnv_df = pd.read_csv(cnv_file, sep='\t')
            
            # Calculate quality metrics
            total_samples = len(cnv_df)
            pass_samples = len(cnv_df[cnv_df['quality_flag'] == 'PASS'])
            warning_samples = len(cnv_df[cnv_df['quality_flag'] == 'WARNING'])
            fail_samples = len(cnv_df[cnv_df['quality_flag'] == 'FAIL'])
            
            avg_quality_score = cnv_df['quality_score'].mean() if 'quality_score' in cnv_df.columns else 0
            
            # Quality thresholds
            pass_rate = pass_samples / total_samples if total_samples > 0 else 0
            overall_status = 'PASS' if pass_rate >= 0.7 and avg_quality_score >= 0.6 else 'WARNING'
            if pass_rate < 0.5 or avg_quality_score < 0.4:
                overall_status = 'FAIL'
            
            self.validation_results['quality_metrics'] = {
                'status': overall_status,
                'total_samples': total_samples,
                'pass_samples': pass_samples,
                'warning_samples': warning_samples,
                'fail_samples': fail_samples,
                'pass_rate': round(pass_rate, 3),
                'average_quality_score': round(avg_quality_score, 3)
            }
            
            return overall_status != 'FAIL'
            
        except Exception as e:
            self.validation_results['quality_metrics'] = {
                'status': 'FAIL',
                'reason': f'Error analyzing quality metrics: {e}'
            }
            return False
    
    def validate_clinical_findings(self):
        """Validate clinical findings make sense"""
        cnv_file = self.results_dir / 'cnv_calls' / 'smn_copy_numbers.txt'
        
        if not cnv_file.exists():
            self.validation_results['clinical_findings'] = {
                'status': 'FAIL',
                'reason': 'CNV calls file not found'
            }
            return False
        
        try:
            cnv_df = pd.read_csv(cnv_file, sep='\t')
            
            # Count clinical findings
            status_counts = cnv_df['sma_status'].value_counts().to_dict()
            
            # Check for logical consistency
            warnings = []
            
            # Check if all samples are AFFECTED (suspicious)
            if status_counts.get('AFFECTED', 0) == len(cnv_df) and len(cnv_df) > 5:
                warnings.append("All samples classified as AFFECTED - check thresholds")
            
            # Check if no carriers found in large cohort (unusual)
            if len(cnv_df) > 50 and status_counts.get('CARRIER', 0) == 0:
                warnings.append("No carriers found in large cohort - unusual pattern")
            
            # Check quality correlation with findings
            if 'quality_score' in cnv_df.columns:
                affected_quality = cnv_df[cnv_df['sma_status'] == 'AFFECTED']['quality_score'].mean()
                if affected_quality < 0.5:
                    warnings.append("Low quality scores for AFFECTED samples - verify findings")
            
            self.validation_results['clinical_findings'] = {
                'status': 'WARNING' if warnings else 'PASS',
                'status_distribution': status_counts,
                'warnings': warnings,
                'total_samples': len(cnv_df)
            }
            
            return True
            
        except Exception as e:
            self.validation_results['clinical_findings'] = {
                'status': 'FAIL',
                'reason': f'Error analyzing clinical findings: {e}'
            }
            return False
    
    def generate_validation_report(self, output_log):
        """Generate comprehensive validation report"""
        with open(output_log, 'w') as f:
            f.write("SMN PIPELINE VALIDATION REPORT\n")
            f.write("=" * 50 + "\n")
            f.write(f"Validation Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Results Directory: {self.results_dir}\n\n")
            
            # Overall status
            all_passed = all(
                result.get('status') in ['PASS', 'WARNING'] 
                for result in self.validation_results.values()
            )
            overall_status = 'PASS' if all_passed else 'FAIL'
            
            f.write(f"OVERALL VALIDATION STATUS: {overall_status}\n\n")
            
            # Detailed results
            for category, results in self.validation_results.items():
                f.write(f"{category.upper().replace('_', ' ')}\n")
                f.write("-" * 30 + "\n")
                f.write(f"Status: {results.get('status', 'UNKNOWN')}\n")
                
                for key, value in results.items():
                    if key != 'status':
                        f.write(f"{key}: {value}\n")
                f.write("\n")
            
            # Recommendations
            f.write("RECOMMENDATIONS\n")
            f.write("-" * 15 + "\n")
            
            if overall_status == 'FAIL':
                f.write("- Pipeline execution issues detected\n")
                f.write("- Review logs for specific errors\n")
                f.write("- Consider re-running failed components\n")
            else:
                f.write("- Pipeline executed successfully\n")
                f.write("- Review individual sample reports\n")
                f.write("- Follow up on critical findings\n")
            
            f.write("\n" + "=" * 50 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Validate SMN pipeline results')
    parser.add_argument('--results-dir', required=True, help='Pipeline results directory')
    parser.add_argument('--expected-samples', type=int, required=True, 
                       help='Expected number of samples')
    parser.add_argument('--critical-exons', required=True,
                       help='Comma-separated list of critical exons')
    parser.add_argument('--output-log', required=True, help='Validation log file')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Initialize validator
    validator = PipelineValidator(args.results_dir)
    
    # Run validations
    logger.info("Starting pipeline validation")
    
    validation_passed = True
    
    # Validate directory structure
    if not validator.validate_directory_structure():
        validation_passed = False
        logger.error("Directory structure validation failed")
    else:
        logger.info("Directory structure validation passed")
    
    # Validate file completeness
    if not validator.validate_file_completeness(args.expected_samples):
        validation_passed = False
        logger.error("File completeness validation failed")
    else:
        logger.info("File completeness validation passed")
    
    # Validate critical exons coverage
    if not validator.validate_critical_exons_coverage(args.critical_exons):
        validation_passed = False
        logger.error("Critical exons validation failed")
    else:
        logger.info("Critical exons validation passed")
    
    # Validate quality metrics
    if not validator.validate_quality_metrics():
        validation_passed = False
        logger.error("Quality metrics validation failed")
    else:
        logger.info("Quality metrics validation passed")
    
    # Validate clinical findings
    if not validator.validate_clinical_findings():
        validation_passed = False
        logger.error("Clinical findings validation failed")
    else:
        logger.info("Clinical findings validation passed")
    
    # Generate validation report
    validator.generate_validation_report(args.output_log)
    
    # Final status
    if validation_passed:
        logger.info("✅ Pipeline validation completed successfully")
        sys.exit(0)
    else:
        logger.error("❌ Pipeline validation failed")
        sys.exit(1)

if __name__ == '__main__':
    main()
