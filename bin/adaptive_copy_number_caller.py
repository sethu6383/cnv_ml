#!/usr/bin/env python3
"""
Adaptive Copy Number Caller for SMN genes
SMA-specific copy number estimation with confidence scoring
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)

class SMACopyNumberCaller:
    """SMA-specific copy number caller with clinical interpretation"""
    
    def __init__(self, threshold_dir):
        self.threshold_dir = Path(threshold_dir)
        self.thresholds = {}
        self.load_thresholds()
        
    def load_thresholds(self):
        """Load current thresholds"""
        threshold_file = self.threshold_dir / "current_thresholds.json"
        
        if threshold_file.exists():
            try:
                with open(threshold_file, 'r') as f:
                    data = json.load(f)
                self.thresholds = data.get('thresholds', {})
                logging.info(f"Loaded thresholds version {data.get('version', 'unknown')}")
            except Exception as e:
                logging.error(f"Error loading thresholds: {e}")
                self._initialize_default_thresholds()
        else:
            self._initialize_default_thresholds()
    
    def _initialize_default_thresholds(self):
        """Initialize default thresholds if none exist"""
        self.thresholds = {
            'SMN1_exon7': {
                'deletion_threshold': -1.5,
                'duplication_threshold': 1.5,
                'copy_number_mapping': {
                    0: {'z_min': -float('inf'), 'z_max': -2.0},
                    1: {'z_min': -2.0, 'z_max': -0.5},
                    2: {'z_min': -0.5, 'z_max': 0.5},
                    3: {'z_min': 0.5, 'z_max': 2.0},
                    4: {'z_min': 2.0, 'z_max': float('inf')}
                }
            },
            'SMN1_exon8': {
                'deletion_threshold': -1.5,
                'duplication_threshold': 1.5,
                'copy_number_mapping': {
                    0: {'z_min': -float('inf'), 'z_max': -2.0},
                    1: {'z_min': -2.0, 'z_max': -0.5},
                    2: {'z_min': -0.5, 'z_max': 0.5},
                    3: {'z_min': 0.5, 'z_max': 2.0},
                    4: {'z_min': 2.0, 'z_max': float('inf')}
                }
            },
            'SMN2_exon7': {
                'deletion_threshold': -1.8,
                'duplication_threshold': 1.8,
                'copy_number_mapping': {
                    0: {'z_min': -float('inf'), 'z_max': -2.5},
                    1: {'z_min': -2.5, 'z_max': -0.8},
                    2: {'z_min': -0.8, 'z_max': 0.8},
                    3: {'z_min': 0.8, 'z_max': 2.5},
                    4: {'z_min': 2.5, 'z_max': float('inf')}
                }
            },
            'SMN2_exon8': {
                'deletion_threshold': -1.8,
                'duplication_threshold': 1.8,
                'copy_number_mapping': {
                    0: {'z_min': -float('inf'), 'z_max': -2.5},
                    1: {'z_min': -2.5, 'z_max': -0.8},
                    2: {'z_min': -0.8, 'z_max': 0.8},
                    3: {'z_min': 0.8, 'z_max': 2.5},
                    4: {'z_min': 2.5, 'z_max': float('inf')}
                }
            }
        }
        logging.info("Initialized default SMA thresholds")
    
    def estimate_copy_number(self, z_score, exon):
        """Estimate copy number from z-score"""
        if exon not in self.thresholds:
            return None, 0.0
        
        mapping = self.thresholds[exon]['copy_number_mapping']
        
        for copy_num, bounds in mapping.items():
            if bounds['z_min'] <= z_score < bounds['z_max']:
                # Calculate confidence based on distance from boundary
                confidence = self._calculate_confidence(z_score, bounds)
                return copy_num, confidence
        
        # Fallback for extreme values
        return 0 if z_score < -3 else 4, 0.5
    
    def _calculate_confidence(self, z_score, bounds):
        """Calculate confidence score based on distance from thresholds"""
        center = (bounds['z_min'] + bounds['z_max']) / 2
        if bounds['z_max'] == float('inf'):
            width = abs(z_score - bounds['z_min'])
            confidence = min(1.0, abs(z_score - bounds['z_min']) / 2.0)
        elif bounds['z_min'] == -float('inf'):
            confidence = min(1.0, abs(z_score - bounds['z_max']) / 2.0)
        else:
            width = bounds['z_max'] - bounds['z_min']
            distance_from_boundary = min(
                abs(z_score - bounds['z_min']),
                abs(z_score - bounds['z_max'])
            )
            confidence = min(1.0, distance_from_boundary / (width / 2))
        
        return max(0.1, confidence)  # Minimum confidence of 0.1
    
    def interpret_sma_status(self, smn1_exon7_cn, smn1_exon8_cn, smn2_exon7_cn, smn2_exon8_cn):
        """Interpret SMA status based on SMN copy numbers"""
        # Focus on SMN1 for primary SMA classification
        smn1_avg_cn = (smn1_exon7_cn + smn1_exon8_cn) / 2
        
        # SMA classification based on SMN1 copy number
        if smn1_avg_cn < 0.5:  # Homozygous deletion
            sma_status = "AFFECTED"
            interpretation = "SMN1 homozygous deletion - likely SMA affected"
            severity_prediction = self._predict_sma_severity(smn2_exon7_cn, smn2_exon8_cn)
        elif smn1_avg_cn < 1.5:  # Heterozygous deletion
            sma_status = "CARRIER"
            interpretation = "SMN1 heterozygous deletion - SMA carrier"
            severity_prediction = "Not applicable (carrier)"
        elif smn1_avg_cn < 2.5:  # Normal copy number
            sma_status = "NORMAL"
            interpretation = "Normal SMN1 copy number"
            severity_prediction = "Not applicable (normal)"
        else:  # Potential duplication
            sma_status = "UNCERTAIN"
            interpretation = "SMN1 copy number above normal - requires confirmation"
            severity_prediction = "Not applicable"
        
        return {
            'sma_status': sma_status,
            'interpretation': interpretation,
            'severity_prediction': severity_prediction,
            'smn1_average_cn': round(smn1_avg_cn, 2),
            'smn2_average_cn': round((smn2_exon7_cn + smn2_exon8_cn) / 2, 2)
        }
    
    def _predict_sma_severity(self, smn2_exon7_cn, smn2_exon8_cn):
        """Predict SMA severity based on SMN2 copy number"""
        smn2_avg_cn = (smn2_exon7_cn + smn2_exon8_cn) / 2
        
        if smn2_avg_cn < 1.5:
            return "Type I (severe) - very low SMN2 copy number"
        elif smn2_avg_cn < 2.5:
            return "Type II-III (intermediate) - low SMN2 copy number"
        elif smn2_avg_cn < 3.5:
            return "Type III-IV (mild) - moderate SMN2 copy number"
        else:
            return "Type IV (mild) - high SMN2 copy number"
    
    def calculate_quality_score(self, confidences):
        """Calculate overall quality score for the sample"""
        if not confidences:
            return 0.0
        
        avg_confidence = np.mean(confidences)
        min_confidence = min(confidences)
        
        # Quality score based on average and minimum confidence
        quality_score = (avg_confidence * 0.7) + (min_confidence * 0.3)
        
        return quality_score
    
    def assign_quality_flag(self, quality_score, exon_consistency):
        """Assign quality flag based on various metrics"""
        if quality_score >= 0.8 and exon_consistency:
            return "PASS"
        elif quality_score >= 0.5:
            return "WARNING"
        else:
            return "FAIL"
    
    def check_exon_consistency(self, smn1_exon7_cn, smn1_exon8_cn, smn2_exon7_cn, smn2_exon8_cn):
        """Check consistency between exon 7 and 8 copy numbers"""
        smn1_diff = abs(smn1_exon7_cn - smn1_exon8_cn)
        smn2_diff = abs(smn2_exon7_cn - smn2_exon8_cn)
        
        # Allow for some variation but flag major discrepancies
        smn1_consistent = smn1_diff <= 1.0
        smn2_consistent = smn2_diff <= 1.5  # SMN2 can be more variable
        
        return smn1_consistent and smn2_consistent

def create_copy_number_plots(df, output_dir):
    """Create visualization plots for copy number results"""
    plt.style.use('default')
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    exons = ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']
    
    for i, exon in enumerate(exons):
        ax = axes[i//2, i%2]
        
        if f"{exon}_z_score" in df.columns and f"{exon}_copy_number" in df.columns:
            # Scatter plot of z-scores vs copy numbers
            scatter = ax.scatter(df[f"{exon}_z_score"], df[f"{exon}_copy_number"], 
                               c=df[f"{exon}_confidence"], cmap='RdYlGn', 
                               alpha=0.7, s=60)
            
            ax.set_xlabel(f"{exon} Z-score")
            ax.set_ylabel(f"{exon} Copy Number")
            ax.set_title(f"{exon} Copy Number Estimation")
            ax.grid(True, alpha=0.3)
            
            # Add colorbar for confidence
            plt.colorbar(scatter, ax=ax, label='Confidence')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "copy_number_scatter.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # SMA status distribution
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    if 'sma_status' in df.columns:
        status_counts = df['sma_status'].value_counts()
        colors = {'NORMAL': 'green', 'CARRIER': 'orange', 'AFFECTED': 'red', 'UNCERTAIN': 'gray'}
        bars = ax.bar(status_counts.index, status_counts.values, 
                     color=[colors.get(status, 'blue') for status in status_counts.index])
        
        ax.set_title("SMA Status Distribution")
        ax.set_ylabel("Number of Samples")
        
        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "sma_status_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Adaptive SMN copy number calling')
    parser.add_argument('--z-scores', required=True, help='Normalized z-scores file')
    parser.add_argument('--threshold-dir', required=True, help='Threshold directory')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--focus-exons', help='Comma-separated list of focus exons')
    parser.add_argument('--generate-plots', action='store_true', help='Generate visualization plots')
    parser.add_argument('--sma-specific', action='store_true', help='Enable SMA-specific interpretation')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load z-scores
    logger.info("Loading normalized z-scores")
    try:
        z_scores_df = pd.read_csv(args.z_scores, sep='\t')
    except Exception as e:
        logger.error(f"Error loading z-scores: {e}")
        sys.exit(1)
    
    logger.info(f"Loaded z-scores for {len(z_scores_df)} samples")
    
    # Initialize copy number caller
    caller = SMACopyNumberCaller(args.threshold_dir)
    
    # Process each sample
    results = []
    exon_columns = [col for col in z_scores_df.columns if col.startswith(('SMN1_', 'SMN2_'))]
    
    for _, row in z_scores_df.iterrows():
        sample_id = row['sample_id']
        sample_type = row.get('sample_type', 'unknown')
        
        logger.debug(f"Processing sample: {sample_id}")
        
        # Estimate copy numbers for each exon
        copy_numbers = {}
        confidences = []
        z_scores = {}
        
        for exon in exon_columns:
            z_score = row[exon]
            copy_num, confidence = caller.estimate_copy_number(z_score, exon)
            
            copy_numbers[exon] = copy_num
            confidences.append(confidence)
            z_scores[exon] = z_score
        
        # SMA-specific interpretation
        if args.sma_specific and all(exon in copy_numbers for exon in 
                                   ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']):
            sma_interpretation = caller.interpret_sma_status(
                copy_numbers['SMN1_exon7'], copy_numbers['SMN1_exon8'],
                copy_numbers['SMN2_exon7'], copy_numbers['SMN2_exon8']
            )
        else:
            sma_interpretation = {
                'sma_status': 'UNKNOWN',
                'interpretation': 'SMA interpretation not performed',
                'severity_prediction': 'Not applicable',
                'smn1_average_cn': 'N/A',
                'smn2_average_cn': 'N/A'
            }
        
        # Quality assessment
        quality_score = caller.calculate_quality_score(confidences)
        exon_consistency = caller.check_exon_consistency(
            copy_numbers.get('SMN1_exon7', 2), copy_numbers.get('SMN1_exon8', 2),
            copy_numbers.get('SMN2_exon7', 2), copy_numbers.get('SMN2_exon8', 2)
        )
        quality_flag = caller.assign_quality_flag(quality_score, exon_consistency)
        
        # Compile result
        result = {
            'sample_id': sample_id,
            'sample_type': sample_type,
            'quality_flag': quality_flag,
            'quality_score': round(quality_score, 3),
            'exon_consistency': exon_consistency,
            **sma_interpretation
        }
        
        # Add individual exon results
        for exon in exon_columns:
            result[f"{exon}_z_score"] = round(z_scores[exon], 3)
            result[f"{exon}_copy_number"] = copy_numbers[exon]
            result[f"{exon}_confidence"] = round(
                confidences[exon_columns.index(exon)], 3
            )
        
        results.append(result)
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results
    output_file = Path(args.output_dir) / "smn_copy_numbers.txt"
    results_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Copy number results saved: {output_file}")
    
    # Generate plots if requested
    if args.generate_plots:
        logger.info("Generating visualization plots")
        create_copy_number_plots(results_df, args.output_dir)
        logger.info("Plots saved to output directory")
    
    # Create summary statistics
    summary = {
        'total_samples': len(results_df),
        'quality_distribution': results_df['quality_flag'].value_counts().to_dict(),
        'sma_status_distribution': results_df['sma_status'].value_counts().to_dict() if 'sma_status' in results_df.columns else {},
        'average_quality_score': float(results_df['quality_score'].mean()),
        'processing_timestamp': pd.Timestamp.now().isoformat()
    }
    
    summary_file = Path(args.output_dir) / "copy_number_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Log summary
    logger.info("Copy number calling completed:")
    logger.info(f"  Total samples: {summary['total_samples']}")
    if 'sma_status_distribution' in summary and summary['sma_status_distribution']:
        logger.info("  SMA status distribution:")
        for status, count in summary['sma_status_distribution'].items():
            logger.info(f"    {status}: {count}")
    logger.info(f"  Average quality score: {summary['average_quality_score']:.3f}")

if __name__ == '__main__':
    main()
            '
