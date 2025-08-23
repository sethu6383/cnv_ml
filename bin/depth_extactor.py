#!/usr/bin/env python3

"""
SMN Depth Extractor with Fallback Analysis
Optimized for SMN1/SMN2 exons 7&8 with intelligent fallback strategies
"""

import os
import sys
import json
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import argparse
import logging
import glob
from dataclasses import dataclass

@dataclass
class ExonDepthResult:
    """Results for depth analysis of a single exon"""
    exon_id: str
    chromosome: str
    start: int
    end: int
    mean_depth: float
    median_depth: float
    positions_covered: int
    total_positions: int
    coverage_fraction: float
    reads_detected: bool
    depth_uniformity: float  # Coefficient of variation
    quality_flags: List[str]

@dataclass
class FallbackRegion:
    """Information about a fallback analysis region"""
    region_id: str
    chromosome: str
    start: int
    end: int
    purpose: str
    read_count: int
    mean_depth: float

@dataclass
class SMNSampleDepth:
    """Complete depth analysis for one sample"""
    sample_id: str
    bam_path: str
    sample_type: str
    analysis_timestamp: str
    primary_exons: Dict[str, ExonDepthResult]
    fallback_regions: List[FallbackRegion]
    fallback_used: bool
    overall_quality: str
    qc_flags: List[str]

class SMNDepthExtractor:
    """Specialized depth extractor for SMN exons with intelligent fallback"""
    
    def __init__(self, bed_file: str, enable_fallback: bool = True):
        self.bed_file = bed_file
        self.enable_fallback = enable_fallback
        
        # Critical exons for SMA
        self.critical_exons = ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']
        
        # Fallback regions for when primary exons have no reads
        self.fallback_regions = {
            'SMN1_upstream': ('chr5', 70945900, 70946065, 'SMN1 upstream region'),
            'SMN1_downstream': ('chr5', 70946177, 70946300, 'SMN1 downstream region'),
            'SMN1_intron7_8': ('chr5', 70946177, 70951940, 'SMN1 intron between exons 7-8'),
            'SMN2_upstream': ('chr5', 70070500, 70070640, 'SMN2 upstream region'),
            'SMN2_downstream': ('chr5', 70070752, 70070900, 'SMN2 downstream region'),
            'SMN2_intron7_8': ('chr5', 70070752, 70076520, 'SMN2 intron between exons 7-8'),
            'SMN_locus_control': ('chr5', 70000000, 70100000, 'SMN locus quality control region'),
        }
        
        # Load exon coordinates
        self.exon_coordinates = self._load_exon_coordinates()
        
    def _load_exon_coordinates(self) -> Dict[str, Tuple[str, int, int]]:
        """Load exon coordinates from BED file"""
        coordinates = {}
        
        try:
            with open(self.bed_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        chrom, start, end, exon_name = parts[:4]
                        coordinates[exon_name] = (chrom, int(start), int(end))
                        
        except Exception as e:
            logging.error(f"Error reading BED file {self.bed_file}: {e}")
            
        return coordinates
    
    def extract_sample_depth(self, bam_path: str, sample_id: str, sample_type: str = "unknown") -> SMNSampleDepth:
        """Extract depth for all SMN exons with fallback analysis"""
        
        logging.info(f"Processing sample: {sample_id}")
        
        # Extract primary exon depths
        primary_results = {}
        for exon_id in self.critical_exons:
            if exon_id in self.exon_coordinates:
                try:
                    depth_result = self._extract_exon_depth(bam_path, exon_id)
                    primary_results[exon_id] = depth_result
                except Exception as e:
                    logging.warning(f"Failed to extract depth for {exon_id}: {e}")
                    # Create failed result
                    coords = self.exon_coordinates[exon_id]
                    primary_results[exon_id] = ExonDepthResult(
                        exon_id=exon_id,
                        chromosome=coords[0],
                        start=coords[1],
                        end=coords[2],
                        mean_depth=0,
                        median_depth=0,
                        positions_covered=0,
                        total_positions=coords[2] - coords[1] + 1,
                        coverage_fraction=0,
                        reads_detected=False,
                        depth_uniformity=0,
                        quality_flags=[f"Extraction failed: {str(e)}"]
                    )
        
        # Assess if fallback is needed
        needs_fallback = self._assess_fallback_need(primary_results)
        fallback_results = []
        
        if needs_fallback and self.enable_fallback:
            logging.info(f"Primary exon coverage insufficient for {sample_id}, performing fallback analysis")
            fallback_results = self._perform_fallback_analysis(bam_path, primary_results)
        
        # Overall quality assessment
        overall_quality, qc_flags = self._assess_overall_quality(primary_results, fallback_results)
        
        return SMNSampleDepth(
            sample_id=sample_id,
            bam_path=bam_path,
            sample_type=sample_type,
            analysis_timestamp=datetime.now().isoformat(),
            primary_exons=primary_results,
            fallback_regions=fallback_results,
            fallback_used=len(fallback_results) > 0,
            overall_quality=overall_quality,
            qc_flags=qc_flags
        )
    
    def _extract_exon_depth(self, bam_path: str, exon_id: str) -> ExonDepthResult:
        """Extract depth for a specific exon"""
        
        if exon_id not in self.exon_coordinates:
            raise ValueError(f"Unknown exon: {exon_id}")
        
        chrom, start, end = self.exon_coordinates[exon_id]
        total_positions = end - start + 1
        
        # Extract depth using samtools depth
        cmd = [
            'samtools', 'depth',
            '-r', f'{chrom}:{start}-{end}',
            '-q', '20',  # Minimum mapping quality
            '-Q', '20',  # Minimum base quality
            '-a',        # Include all positions (even zero depth)
            bam_path
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            if not result.stdout.strip():
                # No output means no reads in region
                return ExonDepthResult(
                    exon_id=exon_id,
                    chromosome=chrom,
                    start=start,
                    end=end,
                    mean_depth=0,
                    median_depth=0,
                    positions_covered=0,
                    total_positions=total_positions,
                    coverage_fraction=0,
                    reads_detected=False,
                    depth_uniformity=0,
                    quality_flags=["NO_READS_DETECTED"]
                )
            
            # Parse depth output
            depths = []
            positions = set()
            
            for line in result.stdout.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    pos = int(parts[1])
                    depth = int(parts[2])
                    depths.append(depth)
                    if depth > 0:
                        positions.add(pos)
            
            if not depths:
                raise ValueError("No depth data parsed")
            
            # Calculate statistics
            mean_depth = np.mean(depths)
            median_depth = np.median(depths)
            coverage_fraction = len(positions) / total_positions
            depth_uniformity = np.std(depths) / mean_depth if mean_depth > 0 else 0
            
            # Quality flags
            quality_flags = []
            if mean_depth < 10:
                quality_flags.append("LOW_DEPTH")
            if coverage_fraction < 0.8:
                quality_flags.append("INCOMPLETE_COVERAGE")
            if depth_uniformity > 1.0:
                quality_flags.append("UNEVEN_COVERAGE")
            
            return ExonDepthResult(
                exon_id=exon_id,
                chromosome=chrom,
                start=start,
                end=end,
                mean_depth=mean_depth,
                median_depth=median_depth,
                positions_covered=len(positions),
                total_positions=total_positions,
                coverage_fraction=coverage_fraction,
                reads_detected=len(positions) > 0,
                depth_uniformity=depth_uniformity,
                quality_flags=quality_flags
            )
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Samtools depth failed for {exon_id}: {e}")
            raise
    
    def _assess_fallback_need(self, primary_results: Dict[str, ExonDepthResult]) -> bool:
        """Assess if fallback analysis is needed"""
        
        # Count exons with no reads or very low coverage
        no_reads_count = 0
        low_coverage_count = 0
        
        for exon_result in primary_results.values():
            if not exon_result.reads_detected:
                no_reads_count += 1
            elif exon_result.mean_depth < 10:
                low_coverage_count += 1
        
        # Fallback if more than 50% of critical exons have issues
        total_critical = len(self.critical_exons)
        problematic_fraction = (no_reads_count + low_coverage_count) / total_critical
        
        return problematic_fraction >= 0.5
    
    def _perform_fallback_analysis(self, bam_path: str, primary_results: Dict[str, ExonDepthResult]) -> List[FallbackRegion]:
        """Perform fallback analysis on alternative regions"""
        
        fallback_results = []
        
        # Identify which gene regions need fallback
        smn1_needs_fallback = any(
            not primary_results[exon].reads_detected 
            for exon in ['SMN1_exon7', 'SMN1_exon8'] 
            if exon in primary_results
        )
        
        smn2_needs_fallback = any(
            not primary_results[exon].reads_detected 
            for exon in ['SMN2_exon7', 'SMN2_exon8'] 
            if exon in primary_results
        )
        
        # Try fallback regions based on needs
        regions_to_try = []
        
        if smn1_needs_fallback:
            regions_to_try.extend(['SMN1_upstream', 'SMN1_downstream', 'SMN1_intron7_8'])
        
        if smn2_needs_fallback:
            regions_to_try.extend(['SMN2_upstream', 'SMN2_downstream', 'SMN2_intron7_8'])
        
        # Always try the general control region
        regions_to_try.append('SMN_locus_control')
        
        for region_id in regions_to_try:
            if region_id in self.fallback_regions:
                try:
                    fallback_result = self._analyze_fallback_region(bam_path, region_id)
                    if fallback_result.read_count > 0:  # Only include regions with reads
                        fallback_results.append(fallback_result)
                except Exception as e:
                    logging.warning(f"Fallback analysis failed for {region_id}: {e}")
        
        return fallback_results
    
    def _analyze_fallback_region(self, bam_path: str, region_id: str) -> FallbackRegion:
        """Analyze a specific fallback region"""
        
        chrom, start, end, purpose = self.fallback_regions[region_id]
        
        # Count reads in region
        read_count_cmd = ['samtools', 'view', '-c', '-q', '20', bam_path, f'{chrom}:{start}-{end}']
        read_result = subprocess.run(read_count_cmd, capture_output=True, text=True, check=True)
        read_count = int(read_result.stdout.strip())
        
        # Get depth if reads are present
        mean_depth = 0
        if read_count > 0:
            depth_cmd = ['samtools', 'depth', '-r', f'{chrom}:{start}-{end}', '-q', '20', '-Q', '20', bam_path]
            depth_result = subprocess.run(depth_cmd, capture_output=True, text=True)
            
            if depth_result.returncode == 0 and depth_result.stdout.strip():
                depths = [int(line.split('\t')[2]) for line in depth_result.stdout.strip().split('\n')]
                if depths:
                    mean_depth = np.mean(depths)
        
        return FallbackRegion(
            region_id=region_id,
            chromosome=chrom,
            start=start,
            end=end,
            purpose=purpose,
            read_count=read_count,
            mean_depth=mean_depth
        )
    
    def _assess_overall_quality(self, primary_results: Dict[str, ExonDepthResult], 
                               fallback_results: List[FallbackRegion]) -> Tuple[str, List[str]]:
        """Assess overall sample quality and generate QC flags"""
        
        qc_flags = []
        
        # Count issues in primary exons
        no_reads_count = sum(1 for result in primary_results.values() if not result.reads_detected)
        low_coverage_count = sum(1 for result in primary_results.values() 
                               if result.reads_detected and result.mean_depth < 10)
        total_exons = len(primary_results)
        
        # Collect specific quality flags
        for exon_result in primary_results.values():
            for flag in exon_result.quality_flags:
                if flag not in qc_flags:
                    qc_flags.append(f"{exon_result.exon_id}: {flag}")
        
        # Add fallback flags
        if fallback_results:
            qc_flags.append(f"FALLBACK_ANALYSIS_USED: {len(fallback_results)} regions")
            for fb in fallback_results:
                if fb.read_count > 0:
                    qc_flags.append(f"FALLBACK_SUCCESS: {fb.region_id} ({fb.read_count} reads)")
                else:
                    qc_flags.append(f"FALLBACK_FAILED: {fb.region_id} (no reads)")
        
        # Determine overall quality
        if no_reads_count == 0 and low_coverage_count == 0:
            overall_quality = "PASS"
        elif no_reads_count <= total_exons * 0.25 and low_coverage_count <= total_exons * 0.25:
            overall_quality = "WARNING"
        elif fallback_results and any(fb.read_count > 0 for fb in fallback_results):
            overall_quality = "WARNING"  # Fallback saved it
        else:
            overall_quality = "FAIL"
        
        return overall_quality, qc_flags


def process_bam_directory(bam_dir: str, bed_file: str, output_dir: str, 
                         sample_type: str = "auto", enable_fallback: bool = True) -> Dict:
    """Process all BAM files in directory"""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Initialize extractor
    extractor = SMNDepthExtractor(bed_file, enable_fallback)
    
    # Find BAM files
    bam_files = glob.glob(os.path.join(bam_dir, "*.bam"))
    
    if not bam_files:
        raise ValueError(f"No BAM files found in {bam_dir}")
    
    logging.info(f"Found {len(bam_files)} BAM files to process")
    
    # Process each BAM file
    results = {}
    summary_stats = {
        'total_samples': len(bam_files),
        'processed_successfully': 0,
        'failed_processing': 0,
        'fallback_used': 0,
        'quality_distribution': {'PASS': 0, 'WARNING': 0, 'FAIL': 0}
    }
    
    for bam_path in bam_files:
        sample_id = os.path.basename(bam_path).replace('.bam', '')
        
        # Auto-detect sample type if needed
        if sample_type == "auto":
            if any(keyword in sample_id.lower() for keyword in ['ref', 'control', 'normal']):
                detected_type = 'reference'
            else:
                detected_type = 'test'
        else:
            detected_type = sample_type
        
        try:
            # Extract depth for this sample
            sample_result = extractor.extract_sample_depth(bam_path, sample_id, detected_type)
            
            # Save individual sample result
            sample_output_file = output_path / f"{sample_id}_depth_results.json"
            with open(sample_output_file, 'w') as f:
                # Convert dataclass to dict for JSON serialization
                sample_dict = {
                    'sample_id': sample_result.sample_id,
                    'bam_path': sample_result.bam_path,
                    'sample_type': sample_result.sample_type,
                    'analysis_timestamp': sample_result.analysis_timestamp,
                    'primary_exons': {
                        exon_id: {
                            'exon_id': exon_result.exon_id,
                            'chromosome': exon_result.chromosome,
                            'start': exon_result.start,
                            'end': exon_result.end,
                            'mean_depth': exon_result.mean_depth,
                            'median_depth': exon_result.median_depth,
                            'positions_covered': exon_result.positions_covered,
                            'total_positions': exon_result.total_positions,
                            'coverage_fraction': exon_result.coverage_fraction,
                            'reads_detected': exon_result.reads_detected,
                            'depth_uniformity': exon_result.depth_uniformity,
                            'quality_flags': exon_result.quality_flags
                        }
                        for exon_id, exon_result in sample_result.primary_exons.items()
                    },
                    'fallback_regions': [
                        {
                            'region_id': fb.region_id,
                            'chromosome': fb.chromosome,
                            'start': fb.start,
                            'end': fb.end,
                            'purpose': fb.purpose,
                            'read_count': fb.read_count,
                            'mean_depth': fb.mean_depth
                        }
                        for fb in sample_result.fallback_regions
                    ],
                    'fallback_used': sample_result.fallback_used,
                    'overall_quality': sample_result.overall_quality,
                    'qc_flags': sample_result.qc_flags
                }
                json.dump(sample_dict, f, indent=2)
            
            results[sample_id] = sample_result
            summary_stats['processed_successfully'] += 1
            
            if sample_result.fallback_used:
                summary_stats['fallback_used'] += 1
            
            summary_stats['quality_distribution'][sample_result.overall_quality] += 1
            
            logging.info(f"Successfully processed {sample_id} (Quality: {sample_result.overall_quality})")
            
        except Exception as e:
            logging.error(f"Failed to process {sample_id}: {e}")
            summary_stats['failed_processing'] += 1
            continue
    
    # Create batch summary
    summary_file = output_path / "batch_depth_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary_stats, f, indent=2)
    
    # Create TSV summary for easy analysis
    tsv_summary = _create_tsv_summary(results, output_path / "depth_summary.tsv")
    
    logging.info(f"Batch processing completed:")
    logging.info(f"  Successfully processed: {summary_stats['processed_successfully']}")
    logging.info(f"  Failed: {summary_stats['failed_processing']}")
    logging.info(f"  Fallback used: {summary_stats['fallback_used']}")
    logging.info(f"  Quality distribution: {summary_stats['quality_distribution']}")
    
    return {
        'results': results,
        'summary_stats': summary_stats,
        'output_directory': str(output_path)
    }


def _create_tsv_summary(results: Dict[str, SMNSampleDepth], output_file: Path) -> str:
    """Create TSV summary of depth results"""
    
    lines = [
        "# SMN Depth Extraction Summary",
        f"# Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"# Total Samples: {len(results)}",
        "#",
        "sample_id\tsample_type\toverall_quality\tfallback_used\t" + 
        "SMN1_exon7_depth\tSMN1_exon7_reads\tSMN1_exon8_depth\tSMN1_exon8_reads\t" +
        "SMN2_exon7_depth\tSMN2_exon7_reads\tSMN2_exon8_depth\tSMN2_exon8_reads\t" +
        "qc_flags"
    ]
    
    for sample_id, sample_result in results.items():
        line_parts = [
            sample_id,
            sample_result.sample_type,
            sample_result.overall_quality,
            str(sample_result.fallback_used)
        ]
        
        # Add depth data for each critical exon
        for exon_id in ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']:
            if exon_id in sample_result.primary_exons:
                exon_result = sample_result.primary_exons[exon_id]
                line_parts.extend([
                    f"{exon_result.mean_depth:.1f}",
                    str(exon_result.reads_detected)
                ])
            else:
                line_parts.extend(["N/A", "False"])
        
        # Add QC flags
        line_parts.append("; ".join(sample_result.qc_flags))
        
        lines.append("\t".join(line_parts))
    
    # Write TSV file
    with open(output_file, 'w') as f:
        f.write("\n".join(lines))
    
    return str(output_file)


def main():
    """Main function for command-line usage"""
    
    parser = argparse.ArgumentParser(
        description='SMN Depth Extractor with Fallback Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic extraction with auto-detection
  %(prog)s --bam-dir /data/bams --bed-file smn_exons.bed --output-dir results/depth
  
  # Reference samples with fallback disabled
  %(prog)s --bam-dir /ref/bams --bed-file smn_exons.bed --output-dir results/depth \\
           --sample-type reference --no-fallback
  
  # Verbose processing with custom quality thresholds
  %(prog)s --bam-dir /data/bams --bed-file smn_exons.bed --output-dir results/depth \\
           --min-depth 15 --min-mapq 30 --verbose
        '''
    )
    
    parser.add_argument('--bam-dir', required=True,
                       help='Directory containing BAM files')
    parser.add_argument('--bed-file', required=True,
                       help='BED file with SMN exon coordinates')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for depth results')
    parser.add_argument('--sample-type', choices=['reference', 'test', 'auto'], 
                       default='auto',
                       help='Sample type classification (default: auto)')
    parser.add_argument('--enable-fallback', action='store_true', default=True,
                       help='Enable fallback analysis for low-coverage samples')
    parser.add_argument('--no-fallback', dest='enable_fallback', action='store_false',
                       help='Disable fallback analysis')
    parser.add_argument('--min-depth', type=int, default=10,
                       help='Minimum depth threshold for quality assessment')
    parser.add_argument('--min-mapq', type=int, default=20,
                       help='Minimum mapping quality for read filtering')
    parser.add_argument('--min-baseq', type=int, default=20,
                       help='Minimum base quality for depth calculation')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    print(f"SMN Depth Extractor with Fallback Analysis")
    print(f"==========================================")
    print(f"BAM Directory: {args.bam_dir}")
    print(f"BED File: {args.bed_file}")
    print(f"Output Directory: {args.output_dir}")
    print(f"Sample Type: {args.sample_type}")
    print(f"Fallback Analysis: {'Enabled' if args.enable_fallback else 'Disabled'}")
    print(f"Quality Thresholds: depth≥{args.min_depth}, mapq≥{args.min_mapq}, baseq≥{args.min_baseq}")
    print("")
    
    try:
        # Process all BAM files
        results = process_bam_directory(
            bam_dir=args.bam_dir,
            bed_file=args.bed_file,
            output_dir=args.output_dir,
            sample_type=args.sample_type,
            enable_fallback=args.enable_fallback
        )
        
        print("Depth extraction completed successfully!")
        print(f"Results saved to: {results['output_directory']}")
        print(f"Samples processed: {results['summary_stats']['processed_successfully']}")
        print(f"Samples with fallback: {results['summary_stats']['fallback_used']}")
        print(f"Quality distribution: {results['summary_stats']['quality_distribution']}")
        
        # Show samples requiring attention
        for sample_id, sample_result in results['results'].items():
            if sample_result.overall_quality == "FAIL":
                print(f"⚠️  ATTENTION: {sample_id} failed quality control")
            elif sample_result.fallback_used:
                print(f"ℹ️  NOTE: {sample_id} used fallback analysis")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
