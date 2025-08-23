#!/usr/bin/env python3
"""
Enhanced SMN Reporter - Generate comprehensive SMA-focused reports
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
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

class SMNReportGenerator:
    """Generate comprehensive SMN/SMA reports"""
    
    def __init__(self, population_cache_dir=None):
        self.population_cache_dir = Path(population_cache_dir) if population_cache_dir else None
        self.population_data = self.load_population_evidence()
        
    def load_population_evidence(self):
        """Load population evidence from cache"""
        if not self.population_cache_dir or not self.population_cache_dir.exists():
            return {}
        
        evidence = {}
        try:
            # Load ClinVar data
            clinvar_file = self.population_cache_dir / "clinvar_SMN1.json"
            if clinvar_file.exists():
                with open(clinvar_file, 'r') as f:
                    evidence['clinvar'] = json.load(f)
            
            # Load population frequencies
            freq_file = self.population_cache_dir / "population_frequencies.json"
            if freq_file.exists():
                with open(freq_file, 'r') as f:
                    evidence['frequencies'] = json.load(f)
                    
        except Exception as e:
            logging.warning(f"Error loading population evidence: {e}")
        
        return evidence
    
    def generate_text_report(self, sample_data, output_file):
        """Generate detailed text report"""
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("SMN COPY NUMBER ANALYSIS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            # Sample information
            f.write("SAMPLE INFORMATION\n")
            f.write("-" * 20 + "\n")
            f.write(f"Sample ID: {sample_data.get('sample_id', 'N/A')}\n")
            f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Sample Type: {sample_data.get('sample_type', 'N/A')}\n")
            f.write(f"Quality Flag: {sample_data.get('quality_flag', 'N/A')}\n")
            f.write(f"Quality Score: {sample_data.get('quality_score', 'N/A')}\n\n")
            
            # SMA interpretation
            f.write("SMA CLINICAL INTERPRETATION\n")
            f.write("-" * 30 + "\n")
            f.write(f"SMA Status: {sample_data.get('sma_status', 'UNKNOWN')}\n")
            f.write(f"Clinical Interpretation: {sample_data.get('interpretation', 'N/A')}\n")
            f.write(f"Predicted Severity: {sample_data.get('severity_prediction', 'N/A')}\n")
            f.write(f"SMN1 Average Copy Number: {sample_data.get('smn1_average_cn', 'N/A')}\n")
            f.write(f"SMN2 Average Copy Number: {sample_data.get('smn2_average_cn', 'N/A')}\n\n")
            
            # Detailed exon results
            f.write("DETAILED EXON ANALYSIS\n")
            f.write("-" * 25 + "\n")
            exons = ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']
            
            for exon in exons:
                z_score = sample_data.get(f"{exon}_z_score", 'N/A')
                copy_num = sample_data.get(f"{exon}_copy_number", 'N/A')
                confidence = sample_data.get(f"{exon}_confidence", 'N/A')
                
                f.write(f"{exon}:\n")
                f.write(f"  Z-score: {z_score}\n")
                f.write(f"  Copy Number: {copy_num}\n")
                f.write(f"  Confidence: {confidence}\n\n")
            
            # Population evidence
            if self.population_data:
                f.write("POPULATION EVIDENCE\n")
                f.write("-" * 20 + "\n")
                
                if 'clinvar' in self.population_data:
                    f.write("ClinVar Evidence:\n")
                    clinvar_data = self.population_data['clinvar']
                    f.write(f"  SMN1 deletions reported: {clinvar_data.get('deletion_reports', 'N/A')}\n")
                    f.write(f"  Clinical significance: {clinvar_data.get('clinical_significance', 'N/A')}\n\n")
                
                if 'frequencies' in self.population_data:
                    f.write("Population Frequencies:\n")
                    freq_data = self.population_data['frequencies']
                    f.write(f"  SMA carrier frequency: {freq_data.get('sma_carrier_freq', 'N/A')}\n")
                    f.write(f"  SMA incidence: {freq_data.get('sma_incidence', 'N/A')}\n\n")
            
            # Clinical recommendations
            f.write("CLINICAL RECOMMENDATIONS\n")
            f.write("-" * 28 + "\n")
            
            sma_status = sample_data.get('sma_status', 'UNKNOWN')
            quality_flag = sample_data.get('quality_flag', 'UNKNOWN')
            
            if sma_status == 'AFFECTED':
                f.write("CRITICAL FINDING:\n")
                f.write("- SMN1 homozygous deletion detected\n")
                f.write("- High probability of SMA\n")
                f.write("- Urgent genetic counseling recommended\n")
                f.write("- Confirmatory testing with MLPA/qPCR advised\n")
                f.write("- Family cascade testing indicated\n")
            elif sma_status == 'CARRIER':
                f.write("IMPORTANT FINDING:\n")
                f.write("- SMN1 heterozygous deletion detected\n")
                f.write("- SMA carrier status\n")
                f.write("- Genetic counseling recommended\n")
                f.write("- Partner testing advised if planning pregnancy\n")
                f.write("- Consider family screening\n")
            elif sma_status == 'NORMAL':
                f.write("NORMAL RESULT:\n")
                f.write("- Standard SMN1 copy numbers detected\n")
                f.write("- No evidence of SMA risk\n")
                f.write("- Routine clinical management\n")
            else:
                f.write("UNCERTAIN RESULT:\n")
                f.write("- Atypical copy number pattern\n")
                f.write("- Confirmatory testing recommended\n")
                f.write("- Clinical correlation advised\n")
            
            if quality_flag == 'FAIL':
                f.write("\nQUALITY CONCERNS:\n")
                f.write("- Poor sample quality detected\n")
                f.write("- Results may be unreliable\n")
                f.write("- Repeat testing strongly recommended\n")
            elif quality_flag == 'WARNING':
                f.write("\nQUALITY NOTE:\n")
                f.write("- Some quality concerns identified\n")
                f.write("- Results should be interpreted with caution\n")
                f.write("- Consider confirmatory testing\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("This report is for research and clinical decision support.\n")
            f.write("Confirmatory testing and genetic counseling are recommended.\n")
            f.write("=" * 80 + "\n")
    
    def generate_html_report(self, sample_data, output_file):
        """Generate interactive HTML report"""
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN Analysis Report - {sample_data.get('sample_id', 'Unknown')}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                   color: white; padding: 20px; border-radius: 10px; margin-bottom: 20px; }}
        .section {{ margin: 20px 0; padding: 15px; border-left: 4px solid #667eea; 
                   background: #f8f9fa; border-radius: 5px; }}
        .critical {{ border-left-color: #dc3545; background: #f8d7da; }}
        .warning {{ border-left-color: #ffc107; background: #fff3cd; }}
        .normal {{ border-left-color: #28a745; background: #d4edda; }}
        .exon-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); 
                     gap: 15px; margin: 15px 0; }}
        .exon-card {{ padding: 15px; background: white; border-radius: 8px; 
                     box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .metric {{ display: flex; justify-content: space-between; margin: 5px 0; }}
        .status-badge {{ padding: 4px 12px; border-radius: 20px; font-weight: bold; 
                        text-align: center; display: inline-block; }}
        .status-affected {{ background: #dc3545; color: white; }}
        .status-carrier {{ background: #ffc107; color: black; }}
        .status-normal {{ background: #28a745; color: white; }}
        .status-uncertain {{ background: #6c757d; color: white; }}
        table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>SMN Copy Number Analysis Report</h1>
        <h2>Sample: {sample_data.get('sample_id', 'Unknown')}</h2>
        <p>Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="section {self._get_status_class(sample_data.get('sma_status', 'UNKNOWN'))}">
        <h3>🧬 SMA Clinical Interpretation</h3>
        <div class="metric">
            <span><strong>SMA Status:</strong></span>
            <span class="status-badge status-{sample_data.get('sma_status', 'unknown').lower()}">
                {sample_data.get('sma_status', 'UNKNOWN')}
            </span>
        </div>
        <div class="metric">
            <span><strong>Clinical Interpretation:</strong></span>
            <span>{sample_data.get('interpretation', 'N/A')}</span>
        </div>
        <div class="metric">
            <span><strong>Predicted Severity:</strong></span>
            <span>{sample_data.get('severity_prediction', 'N/A')}</span>
        </div>
        <div class="metric">
            <span><strong>SMN1 Average Copy Number:</strong></span>
            <span>{sample_data.get('smn1_average_cn', 'N/A')}</span>
        </div>
        <div class="metric">
            <span><strong>SMN2 Average Copy Number:</strong></span>
            <span>{sample_data.get('smn2_average_cn', 'N/A')}</span>
        </div>
    </div>
    
    <div class="section">
        <h3>📊 Detailed Exon Analysis</h3>
        <div class="exon-grid">
            {self._generate_exon_cards(sample_data)}
        </div>
    </div>
    
    <div class="section">
        <h3>🔬 Quality Metrics</h3>
        <div class="metric">
            <span><strong>Quality Flag:</strong></span>
            <span class="status-badge {self._get_quality_class(sample_data.get('quality_flag', 'UNKNOWN'))}">
                {sample_data.get('quality_flag', 'UNKNOWN')}
            </span>
        </div>
        <div class="metric">
            <span><strong>Quality Score:</strong></span>
            <span>{sample_data.get('quality_score', 'N/A')}</span>
        </div>
        <div class="metric">
            <span><strong>Exon Consistency:</strong></span>
            <span>{'✓ Consistent' if sample_data.get('exon_consistency', False) else '⚠ Inconsistent'}</span>
        </div>
    </div>
    
    {self._generate_population_evidence_section()}
    
    {self._generate_recommendations_section(sample_data)}
    
    <div class="section">
        <h3>📋 Technical Details</h3>
        <table>
            <tr><th>Parameter</th><th>Value</th></tr>
            <tr><td>Sample Type</td><td>{sample_data.get('sample_type', 'N/A')}</td></tr>
            <tr><td>Processing Pipeline</td><td>SMN CNV v2.0-SMA-OPTIMIZED</td></tr>
            <tr><td>Reference Genome</td><td>GRCh38</td></tr>
            <tr><td>Analysis Focus</td><td>SMN1/SMN2 Exons 7 & 8</td></tr>
        </table>
    </div>
    
    <footer style="margin-top: 30px; padding: 20px; background: #f8f9fa; border-radius: 5px;">
        <p><strong>Disclaimer:</strong> This analysis is for research and clinical decision support. 
        Confirmatory testing with validated clinical methods is recommended. 
        Genetic counseling should be provided for actionable findings.</p>
    </footer>
</body>
</html>
"""
        
        with open(output_file, 'w') as f:
            f.write(html_content)
    
    def _get_status_class(self, status):
        """Get CSS class for SMA status"""
        status_classes = {
            'AFFECTED': 'critical',
            'CARRIER': 'warning', 
            'NORMAL': 'normal',
            'UNCERTAIN': 'warning'
        }
        return status_classes.get(status, '')
    
    def _get_quality_class(self, quality):
        """Get CSS class for quality flag"""
        quality_classes = {
            'PASS': 'status-normal',
            'WARNING': 'status-carrier',
            'FAIL': 'status-affected'
        }
        return quality_classes.get(quality, 'status-uncertain')
    
    def _generate_exon_cards(self, sample_data):
        """Generate HTML cards for each exon"""
        exons = ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']
        cards_html = ""
        
        for exon in exons:
            z_score = sample_data.get(f"{exon}_z_score", 'N/A')
            copy_num = sample_data.get(f"{exon}_copy_number", 'N/A')
            confidence = sample_data.get(f"{exon}_confidence", 'N/A')
            
            cards_html += f"""
            <div class="exon-card">
                <h4>{exon}</h4>
                <div class="metric">
                    <span><strong>Z-score:</strong></span>
                    <span>{z_score}</span>
                </div>
                <div class="metric">
                    <span><strong>Copy Number:</strong></span>
                    <span>{copy_num}</span>
                </div>
                <div class="metric">
                    <span><strong>Confidence:</strong></span>
                    <span>{confidence}</span>
                </div>
            </div>
            """
        
        return cards_html
    
    def _generate_population_evidence_section(self):
        """Generate population evidence section"""
        if not self.population_data:
            return ""
        
        html = """
    <div class="section">
        <h3>🌍 Population Evidence</h3>
        """
        
        if 'clinvar' in self.population_data:
            html += """
        <h4>ClinVar Database</h4>
        <p>SMN1 deletions are well-documented in ClinVar with strong clinical evidence for SMA causation.</p>
        """
        
        if 'frequencies' in self.population_data:
            html += """
        <h4>Population Frequencies</h4>
        <p>SMA carrier frequency: ~1 in 50 individuals globally<br>
        SMA incidence: ~1 in 10,000 live births</p>
        """
        
        html += "</div>"
        return html
    
    def _generate_recommendations_section(self, sample_data):
        """Generate clinical recommendations section"""
        sma_status = sample_data.get('sma_status', 'UNKNOWN')
        quality_flag = sample_data.get('quality_flag', 'UNKNOWN')
        
        if sma_status == 'AFFECTED':
            rec_class = 'critical'
            rec_title = '🚨 Critical Clinical Action Required'
            recommendations = [
                "Immediate genetic counseling",
                "Confirmatory testing with MLPA or qPCR",
                "Clinical correlation with symptoms",
                "Family cascade testing",
                "Consider pediatric neurology referral"
            ]
        elif sma_status == 'CARRIER':
            rec_class = 'warning'
            rec_title = '⚠️ Important Clinical Considerations'
            recommendations = [
                "Genetic counseling recommended",
                "Partner screening if planning pregnancy",
                "Confirmatory testing advised",
                "Family screening considerations",
                "Reproductive planning counseling"
            ]
        elif sma_status == 'NORMAL':
            rec_class = 'normal'
            rec_title = '✅ Normal Result - Routine Management'
            recommendations = [
                "No immediate clinical action required",
                "Routine clinical management",
                "Consider in differential diagnosis if symptoms present"
            ]
        else:
            rec_class = 'warning'
            rec_title = '❓ Uncertain Result - Further Investigation'
            recommendations = [
                "Confirmatory testing strongly recommended",
                "Clinical correlation essential",
                "Consider repeat sampling",
                "Genetic counseling may be helpful"
            ]
        
        if quality_flag in ['WARNING', 'FAIL']:
            recommendations.insert(0, "Address quality concerns with repeat testing")
        
        rec_html = f"""
    <div class="section {rec_class}">
        <h3>{rec_title}</h3>
        <ul>
        """
        
        for rec in recommendations:
            rec_html += f"<li>{rec}</li>"
        
        rec_html += """
        </ul>
    </div>
        """
        
        return rec_html

def load_cnv_results(cnv_file):
    """Load CNV calling results"""
    try:
        return pd.read_csv(cnv_file, sep='\t')
    except Exception as e:
        logging.error(f"Error loading CNV results: {e}")
        return pd.DataFrame()

def load_depth_results(depth_dir):
    """Load depth extraction results"""
    depth_dir = Path(depth_dir)
    depth_data = {}
    
    for json_file in depth_dir.glob("*_depth_results.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            depth_data[data['sample_id']] = data
        except Exception as e:
            logging.warning(f"Error loading {json_file}: {e}")
    
    return depth_data

def main():
    parser = argparse.ArgumentParser(description='Generate enhanced SMN reports')
    parser.add_argument('--cnv-results', required=True, help='CNV calling results file')
    parser.add_argument('--depth-results', required=True, help='Depth results directory')
    parser.add_argument('--population-cache', help='Population evidence cache directory')
    parser.add_argument('--threshold-history', help='Threshold history directory')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--format', default='all', choices=['txt', 'html', 'json', 'all'],
                       help='Report format')
    parser.add_argument('--include-population-evidence', action='store_true',
                       help='Include population evidence in reports')
    parser.add_argument('--sma-clinical-context', action='store_true',
                       help='Include SMA clinical context')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load CNV results
    logger.info("Loading CNV calling results")
    cnv_df = load_cnv_results(args.cnv_results)
    
    if cnv_df.empty:
        logger.error("No CNV results found")
        sys.exit(1)
    
    # Load depth results
    logger.info("Loading depth extraction results")
    depth_data = load_depth_results(args.depth_results)
    
    # Initialize report generator
    reporter = SMNReportGenerator(args.population_cache)
    
    # Generate reports for each sample
    logger.info(f"Generating reports for {len(cnv_df)} samples")
    
    for _, row in cnv_df.iterrows():
        sample_id = row['sample_id']
        sample_dir = Path(args.output_dir) / sample_id
        sample_dir.mkdir(exist_ok=True)
        
        # Combine CNV and depth data
        sample_data = row.to_dict()
        if sample_id in depth_data:
            sample_data.update(depth_data[sample_id])
        
        # Generate requested formats
        if args.format in ['txt', 'all']:
            txt_file = sample_dir / f"{sample_id}_report.txt"
            reporter.generate_text_report(sample_data, txt_file)
        
        if args.format in ['html', 'all']:
            html_file = sample_dir / f"{sample_id}_report.html"
            reporter.generate_html_report(sample_data, html_file)
        
        if args.format in ['json', 'all']:
            json_file = sample_dir / f"{sample_id}_data.json"
            with open(json_file, 'w') as f:
                json.dump(sample_data, f, indent=2, default=str)
    
    logger.info("Enhanced SMN report generation completed")

if __name__ == '__main__':
    main()
