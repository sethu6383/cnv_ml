#!/usr/bin/env python3
"""
SMN Depth Extractor - Extract read depth for SMN critical exons
Focus on SMN1/SMN2 exons 7&8 with fallback analysis
"""

import os
import sys
import json
import argparse
import subprocess
import pandas as pd
from pathlib import Path
import logging
from collections import defaultdict

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)

def parse_bed_file(bed_file):
    """Parse BED file to extract target regions"""
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    chrom, start, end, name = parts[:4]
                    regions.append({
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'name': name
                    })
    return regions

def extract_depth_samtools(bam_file, regions, min_mapq=20, min_baseq=20):
    """Extract depth using samtools for each region"""
    depth_data = {}
    
    for region in regions:
        region_str = f"{region['chrom']}:{region['start']}-{region['end']}"
        
        cmd = [
            'samtools', 'depth',
            '-r', region_str,
            '-q', str(min_mapq),
            '-Q', str(min_baseq),
            bam_file
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            depths = []
            
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        depths.append(int(parts[2]))
            
            if depths:
                depth_data[region['name']] = {
                    'mean_depth': sum(depths) / len(depths),
                    'median_depth': sorted(depths)[len(depths)//2],
                    'min_depth': min(depths),
                    'max_depth': max(depths),
                    'coverage_breadth': len(depths),
                    'total_positions': region['end'] - region['start'],
                    'coverage_fraction': len(depths) / (region['end'] - region['start'])
                }
            else:
                depth_data[region['name']] = {
                    'mean_depth': 0,
                    'median_depth': 0,
                    'min_depth': 0,
                    'max_depth': 0,
                    'coverage_breadth': 0,
                    'total_positions': region['end'] - region['start'],
                    'coverage_fraction': 0
                }
                
        except subprocess.CalledProcessError as e:
            logging.error(f"Error extracting depth for {region['name']}: {e}")
            depth_data[region['name']] = None
    
    return depth_data

def assess_sample_quality(depth_data, min_depth=10):
    """Assess overall sample quality based on depth metrics"""
    total_regions = len(depth_data)
    passed_regions = 0
    total_coverage = 0
    
    for region_name, data in depth_data.items():
        if data and data['mean_depth'] >= min_depth and data['coverage_fraction'] >= 0.8:
            passed_regions += 1
        if data:
            total_coverage += data['coverage_fraction']
    
    quality_score = passed_regions / total_regions if total_regions > 0 else 0
    avg_coverage_fraction = total_coverage / total_regions if total_regions > 0 else 0
    
    if quality_score >= 0.8 and avg_coverage_fraction >= 0.8:
        quality = "PASS"
    elif quality_score >= 0.5:
        quality = "WARNING"
    else:
        quality = "FAIL"
    
    return {
        'quality': quality,
        'quality_score': quality_score,
        'regions_passed': passed_regions,
        'total_regions': total_regions,
        'avg_coverage_fraction': avg_coverage_fraction
    }

def determine_sample_type(filename, sample_type_arg):
    """Determine sample type from filename or argument"""
    if sample_type_arg != 'auto':
        return sample_type_arg
    
    filename_lower = filename.lower()
    if any(keyword in filename_lower for keyword in ['ref', 'control', 'normal']):
        return 'reference'
    else:
        return 'test'

def process_bam_file(bam_file, regions, output_dir, sample_type_arg, min_depth, min_mapq, min_baseq):
    """Process a single BAM file"""
    logger = logging.getLogger(__name__)
    sample_id = Path(bam_file).stem.replace('.bam', '')
    
    logger.info(f"Processing sample: {sample_id}")
    
    # Check if BAM file exists and is indexed
    if not os.path.exists(bam_file):
        logger.error(f"BAM file not found: {bam_file}")
        return None
    
    index_file = f"{bam_file}.bai"
    if not os.path.exists(index_file):
        logger.warning(f"BAM index not found, creating: {index_file}")
        try:
            subprocess.run(['samtools', 'index', bam_file], check=True)
        except subprocess.CalledProcessError:
            logger.error(f"Failed to create index for {bam_file}")
            return None
    
    # Extract depth data
    depth_data = extract_depth_samtools(bam_file, regions, min_mapq, min_baseq)
    
    # Assess quality
    quality_metrics = assess_sample_quality(depth_data, min_depth)
    
    # Determine sample type
    sample_type = determine_sample_type(sample_id, sample_type_arg)
    
    # Compile results
    results = {
        'sample_id': sample_id,
        'bam_file': bam_file,
        'sample_type': sample_type,
        'quality_metrics': quality_metrics,
        'depth_data': depth_data,
        'processing_parameters': {
            'min_depth': min_depth,
            'min_mapq': min_mapq,
            'min_baseq': min_baseq,
            'regions_analyzed': len(regions)
        },
        'timestamp': pd.Timestamp.now().isoformat()
    }
    
    # Save results
    output_file = Path(output_dir) / f"{sample_id}_depth_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Results saved: {output_file}")
    return results

def main():
    parser = argparse.ArgumentParser(description='Extract SMN gene depth from BAM files')
    parser.add_argument('--bam-dir', required=True, help='Directory containing BAM files')
    parser.add_argument('--bed-file', required=True, help='BED file with target regions')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--sample-type', default='auto', choices=['reference', 'test', 'auto'],
                       help='Sample type classification')
    parser.add_argument('--enable-fallback', action='store_true',
                       help='Enable fallback analysis for low coverage')
    parser.add_argument('--min-depth', type=int, default=10,
                       help='Minimum depth threshold')
    parser.add_argument('--min-mapq', type=int, default=20,
                       help='Minimum mapping quality')
    parser.add_argument('--min-baseq', type=int, default=20,
                       help='Minimum base quality')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse BED file
    logger.info(f"Parsing BED file: {args.bed_file}")
    regions = parse_bed_file(args.bed_file)
    logger.info(f"Found {len(regions)} target regions")
    
    # Find BAM files
    bam_files = list(Path(args.bam_dir).glob('*.bam'))
    if not bam_files:
        logger.error(f"No BAM files found in {args.bam_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(bam_files)} BAM files")
    
    # Process each BAM file
    results_summary = []
    for bam_file in bam_files:
        try:
            result = process_bam_file(
                str(bam_file), regions, args.output_dir,
                args.sample_type, args.min_depth, args.min_mapq, args.min_baseq
            )
            if result:
                results_summary.append({
                    'sample_id': result['sample_id'],
                    'sample_type': result['sample_type'],
                    'quality': result['quality_metrics']['quality'],
                    'quality_score': result['quality_metrics']['quality_score'],
                    'regions_passed': result['quality_metrics']['regions_passed']
                })
        except Exception as e:
            logger.error(f"Error processing {bam_file}: {e}")
    
    # Save summary
    if results_summary:
        summary_df = pd.DataFrame(results_summary)
        summary_file = Path(args.output_dir) / 'depth_extraction_summary.tsv'
        summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"Summary saved: {summary_file}")
        
        # Log summary statistics
        logger.info(f"Processing complete:")
        logger.info(f"  Total samples: {len(results_summary)}")
        logger.info(f"  Quality PASS: {len(summary_df[summary_df['quality'] == 'PASS'])}")
        logger.info(f"  Quality WARNING: {len(summary_df[summary_df['quality'] == 'WARNING'])}")
        logger.info(f"  Quality FAIL: {len(summary_df[summary_df['quality'] == 'FAIL'])}")
    
    logger.info("SMN depth extraction completed")

if __name__ == '__main__':
    main()
