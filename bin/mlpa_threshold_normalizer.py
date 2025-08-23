#!/usr/bin/env python3
"""
MLPA Threshold Normalizer - Adaptive threshold learning from MLPA data
Normalize SMN depth data with population-validated thresholds
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import matthews_corrcoef
import pickle

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)

class MLPAThresholdManager:
    """Manages adaptive thresholds learned from MLPA training data"""
    
    def __init__(self, threshold_dir, mlpa_file=None):
        self.threshold_dir = Path(threshold_dir)
        self.threshold_dir.mkdir(exist_ok=True)
        self.mlpa_file = mlpa_file
        self.thresholds = {}
        self.performance_history = []
        self.version = "1.0"
        
    def load_mlpa_training_data(self):
        """Load MLPA training dataset"""
        if not self.mlpa_file or not os.path.exists(self.mlpa_file):
            logging.warning("No MLPA training file provided or found")
            return None
        
        try:
            # Expected format: sample_id, SMN1_copy_number, SMN2_copy_number, clinical_interpretation
            df = pd.read_csv(self.mlpa_file, sep='\t')
            logging.info(f"Loaded {len(df)} MLPA training samples")
            return df
        except Exception as e:
            logging.error(f"Error loading MLPA data: {e}")
            return None
    
    def load_existing_thresholds(self):
        """Load existing threshold configuration"""
        threshold_file = self.threshold_dir / "current_thresholds.json"
        
        if threshold_file.exists():
            try:
                with open(threshold_file, 'r') as f:
                    data = json.load(f)
                self.thresholds = data.get('thresholds', {})
                self.version = data.get('version', '1.0')
                self.performance_history = data.get('performance_history', [])
                logging.info(f"Loaded existing thresholds version {self.version}")
                return True
            except Exception as e:
                logging.error(f"Error loading thresholds: {e}")
        
        return False
    
    def initialize_default_thresholds(self):
        """Initialize default SMA-optimized thresholds"""
        self.thresholds = {
            'SMN1_exon7': {
                'deletion_threshold': -1.5,  # Z-score threshold for deletion
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
                'deletion_threshold': -1.8,  # SMN2 typically more variable
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
        
        logging.info("Initialized default SMA-optimized thresholds")
    
    def optimize_thresholds_with_mlpa(self, depth_data, mlpa_data):
        """Optimize thresholds using MLPA training data"""
        if mlpa_data is None or depth_data.empty:
            logging.warning("Insufficient data for threshold optimization")
            return
        
        # Match samples between depth data and MLPA data
        common_samples = set(depth_data['sample_id']) & set(mlpa_data['sample_id'])
        if len(common_samples) < 10:
            logging.warning(f"Too few common samples ({len(common_samples)}) for threshold optimization")
            return
        
        logging.info(f"Optimizing thresholds with {len(common_samples)} training samples")
        
        # For each exon, optimize thresholds
        for exon in ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']:
            if exon not in depth_data.columns:
                continue
                
            self._optimize_exon_thresholds(depth_data, mlpa_data, exon, common_samples)
        
        # Update version
        self.version = f"{float(self.version) + 0.1:.1f}"
        logging.info(f"Thresholds optimized, new version: {self.version}")
    
    def _optimize_exon_thresholds(self, depth_data, mlpa_data, exon, common_samples):
        """Optimize thresholds for a specific exon"""
        # Extract relevant data
        depth_subset = depth_data[depth_data['sample_id'].isin(common_samples)]
        mlpa_subset = mlpa_data[mlpa_data['sample_id'].isin(common_samples)]
        
        # Get expected copy numbers from MLPA
        gene_name = exon.split('_')[0]  # SMN1 or SMN2
        mlpa_col = f"{gene_name}_copy_number"
        
        if mlpa_col not in mlpa_subset.columns:
            logging.warning(f"MLPA column {mlpa_col} not found")
            return
        
        # Merge data
        merged = pd.merge(depth_subset[['sample_id', exon]], 
                         mlpa_subset[['sample_id', mlpa_col]], 
                         on='sample_id')
        
        if len(merged) < 5:
            return
        
        # Calculate optimal thresholds using ROC analysis
        z_scores = merged[exon].values
        true_copy_numbers = merged[mlpa_col].values
        
        # Optimize deletion threshold (0 vs 1+ copies)
        deletion_labels = (true_copy_numbers == 0).astype(int)
        if sum(deletion_labels) > 0 and sum(deletion_labels) < len(deletion_labels):
            del_threshold = self._find_optimal_threshold(z_scores, deletion_labels, mode='deletion')
            if del_threshold is not None:
                self.thresholds[exon]['deletion_threshold'] = del_threshold
        
        # Optimize duplication threshold (3+ vs 2 copies)
        dup_labels = (true_copy_numbers >= 3).astype(int)
        if sum(dup_labels) > 0 and sum(dup_labels) < len(dup_labels):
            dup_threshold = self._find_optimal_threshold(z_scores, dup_labels, mode='duplication')
            if dup_threshold is not None:
                self.thresholds[exon]['duplication_threshold'] = dup_threshold
    
    def _find_optimal_threshold(self, scores, labels, mode='deletion'):
        """Find optimal threshold using Matthews Correlation Coefficient"""
        best_threshold = None
        best_mcc = -1
        
        # Test thresholds in reasonable range
        if mode == 'deletion':
            test_thresholds = np.arange(-3.0, 0.0, 0.1)
        else:  # duplication
            test_thresholds = np.arange(0.0, 3.0, 0.1)
        
        for threshold in test_thresholds:
            if mode == 'deletion':
                predictions = (scores < threshold).astype(int)
            else:
                predictions = (scores > threshold).astype(int)
            
            mcc = matthews_corrcoef(labels, predictions)
            if mcc > best_mcc:
                best_mcc = mcc
                best_threshold = threshold
        
        return best_threshold if best_mcc > 0.1 else None
    
    def save_thresholds(self):
        """Save current thresholds to file"""
        threshold_data = {
            'version': self.version,
            'thresholds': self.thresholds,
            'performance_history': self.performance_history,
            'timestamp': pd.Timestamp.now().isoformat()
        }
        
        # Save current version
        current_file = self.threshold_dir / "current_thresholds.json"
        with open(current_file, 'w') as f:
            json.dump(threshold_data, f, indent=2)
        
        # Save versioned backup
        version_file = self.threshold_dir / f"thresholds_v{self.version}.json"
        with open(version_file, 'w') as f:
            json.dump(threshold_data, f, indent=2)
        
        logging.info(f"Thresholds saved: version {self.version}")

def load_depth_data(depth_dir):
    """Load all depth result files"""
    depth_dir = Path(depth_dir)
    all_data = []
    
    for json_file in depth_dir.glob("*_depth_results.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            sample_id = data['sample_id']
            sample_type = data['sample_type']
            
            # Extract depth metrics for each exon
            row = {'sample_id': sample_id, 'sample_type': sample_type}
            
            for exon_name, exon_data in data['depth_data'].items():
                if exon_data:
                    row[exon_name] = exon_data['mean_depth']
                else:
                    row[exon_name] = 0
            
            all_data.append(row)
            
        except Exception as e:
            logging.error(f"Error loading {json_file}: {e}")
    
    return pd.DataFrame(all_data)

def normalize_depth_data(depth_df):
    """Normalize depth data using reference samples"""
    # Separate reference and test samples
    ref_samples = depth_df[depth_df['sample_type'] == 'reference']
    
    if len(ref_samples) < 3:
        logging.warning("Insufficient reference samples, using all samples for normalization")
        ref_samples = depth_df
    
    # Calculate z-scores for each exon
    exon_columns = [col for col in depth_df.columns if col.startswith(('SMN1_', 'SMN2_'))]
    normalized_df = depth_df.copy()
    
    for exon in exon_columns:
        ref_mean = ref_samples[exon].mean()
        ref_std = ref_samples[exon].std()
        
        if ref_std > 0:
            normalized_df[exon] = (depth_df[exon] - ref_mean) / ref_std
        else:
            normalized_df[exon] = 0
    
    return normalized_df

def main():
    parser = argparse.ArgumentParser(description='MLPA-optimized threshold normalization')
    parser.add_argument('--depth-dir', required=True, help='Directory with depth results')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--threshold-dir', required=True, help='Threshold management directory')
    parser.add_argument('--mlpa-file', help='MLPA training data file')
    parser.add_argument('--population-cache', help='Population evidence cache directory')
    parser.add_argument('--force-retrain', action='store_true', help='Force threshold retraining')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize threshold manager
    threshold_manager = MLPAThresholdManager(args.threshold_dir, args.mlpa_file)
    
    # Load or initialize thresholds
    if not args.force_retrain and threshold_manager.load_existing_thresholds():
        logger.info("Using existing thresholds")
    else:
        threshold_manager.initialize_default_thresholds()
    
    # Load depth data
    logger.info("Loading depth extraction results")
    depth_df = load_depth_data(args.depth_dir)
    
    if depth_df.empty:
        logger.error("No depth data found")
        sys.exit(1)
    
    logger.info(f"Loaded depth data for {len(depth_df)} samples")
    
    # Load MLPA training data if available
    mlpa_data = threshold_manager.load_mlpa_training_data()
    
    # Optimize thresholds with MLPA data
    if mlpa_data is not None and (args.force_retrain or len(threshold_manager.performance_history) == 0):
        threshold_manager.optimize_thresholds_with_mlpa(depth_df, mlpa_data)
    
    # Normalize depth data
    logger.info("Normalizing depth data")
    normalized_df = normalize_depth_data(depth_df)
    
    # Save normalized data
    output_file = Path(args.output_dir) / "z_scores_optimized.txt"
    normalized_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Normalized data saved: {output_file}")
    
    # Save thresholds
    threshold_manager.save_thresholds()
    
    # Create summary
    exon_columns = [col for col in normalized_df.columns if col.startswith(('SMN1_', 'SMN2_'))]
    summary = {
        'total_samples': len(normalized_df),
        'reference_samples': len(normalized_df[normalized_df['sample_type'] == 'reference']),
        'test_samples': len(normalized_df[normalized_df['sample_type'] == 'test']),
        'exons_analyzed': exon_columns,
        'threshold_version': threshold_manager.version,
        'normalization_stats': {}
    }
    
    for exon in exon_columns:
        summary['normalization_stats'][exon] = {
            'mean': float(normalized_df[exon].mean()),
            'std': float(normalized_df[exon].std()),
            'min': float(normalized_df[exon].min()),
            'max': float(normalized_df[exon].max())
        }
    
    summary_file = Path(args.output_dir) / "normalization_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info("MLPA threshold normalization completed")

if __name__ == '__main__':
    main()
