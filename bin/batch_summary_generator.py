#!/usr/bin/env python3
"""
Batch Summary Generator - Create comprehensive batch analysis summaries
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
import matplotlib.pyplot as plt
import seaborn as sns

def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)

class BatchSummaryGenerator:
    """Generate batch analysis summaries and dashboards"""
    
    def __init__(self):
        self.timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
    def load_sample_reports(self, reports_dir):
        """Load all sample report data"""
        reports_dir = Path(reports_dir)
        all_data = []
        
        for sample_dir in reports_dir.iterdir():
            if sample_dir.is_dir():
                json_file = sample_dir / f"{sample_dir.name}_data.json"
                if json_file.exists():
                    try:
                        with open(json_file, 'r') as f:
                            data = json.load(f)
                        all_data.append(data)
                    except Exception as e:
                        logging.warning(f"Error loading {json_file}: {e}")
        
        return pd.DataFrame(all_data)
    
    def generate_batch_tsv(self, df, output_file):
        """Generate tabular batch summary"""
        # Select key columns for summary
        summary_columns = [
            'sample_id', 'sample_type', 'sma_status', 'quality_flag',
            'quality_score', 'smn1_average_cn', 'smn2_average_cn',
            'interpretation', 'severity_prediction', 'exon_consistency'
        ]
        
        # Add exon-specific data
        exon_columns = []
        for exon in ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']:
            for metric in ['z_score', 'copy_number', 'confidence']:
                col_name = f"{exon}_{metric}"
                if col_name in df.columns:
                    summary_columns.append(col_name)
                    exon_columns.append(col_name)
        
        # Create summary DataFrame
        summary_df = df[summary_columns].copy()
        
        # Add derived metrics
        summary_df['critical_finding'] = summary_df['sma_status'].isin(['AFFECTED', 'CARRIER'])
        summary_df['needs_followup'] = (
            summary_df['sma_status'].isin(['AFFECTED', 'CARRIER', 'UNCERTAIN']) |
            summary_df['quality_flag'].isin(['WARNING', 'FAIL'])
        )
        
        # Save summary
        summary_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Batch TSV summary saved: {output_file}")
        
        return summary_df
    
    def generate_batch_html_dashboard(self, df, output_file):
        """Generate interactive HTML dashboard"""
        # Calculate summary statistics
        total_samples = len(df)
        status_counts = df['sma_status'].value_counts().to_dict()
        quality_counts = df['quality_flag'].value_counts().to_dict()
        
        # Generate statistics table
        stats_html = self._generate_statistics_table(df)
        
        # Generate sample table
        table_html = self._generate_sample_table(df)
        
        # Generate plots as base64 if matplotlib available
        plots_html = self._generate_plots_html(df)
        
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN Batch Analysis Dashboard</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
               margin: 0; padding: 20px; background: #f5f5f5; }}
        .dashboard {{ max-width: 1200px; margin: 0 auto; }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                   color: white; padding: 30px; border-radius: 15px; margin-bottom: 30px; text-align: center; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); 
                      gap: 20px; margin: 30px 0; }}
        .stat-card {{ background: white; padding: 20px; border-radius: 10px; 
                     box-shadow: 0 4px 6px rgba(0,0,0,0.1); text-align: center; }}
        .stat-number {{ font-size: 2.5em; font-weight: bold; margin: 10px 0; }}
        .stat-label {{ color: #666; font-size: 0.9em; }}
        .critical {{ color: #dc3545; }}
        .warning {{ color: #ffc107; }}
        .success {{ color: #28a745; }}
        .info {{ color: #17a2b8; }}
        .section {{ background: white; margin: 20px 0; padding: 25px; 
                   border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background: #f8f9fa; font-weight: 600; }}
        .status-affected {{ background: #dc3545; color: white; padding: 4px 8px; border-radius: 4px; }}
        .status-carrier {{ background: #ffc107; color: black; padding: 4px 8px; border-radius: 4px; }}
        .status-normal {{ background: #28a745; color: white; padding: 4px 8px; border-radius: 4px; }}
        .status-uncertain {{ background: #6c757d; color: white; padding: 4px 8px; border-radius: 4px; }}
        .quality-pass {{ background: #28a745; color: white; padding: 2px 6px; border-radius: 3px; }}
        .quality-warning {{ background: #ffc107; color: black; padding: 2px 6px; border-radius: 3px; }}
        .quality-fail {{ background: #dc3545; color: white; padding: 2px 6px; border-radius: 3px; }}
        .plots {{ margin: 20px 0; text-align: center; }}
        .footer {{ text-align: center; margin-top: 40px; padding: 20px; 
                  background: #e9ecef; border-radius: 10px; }}
    </style>
    <script>
        function filterTable() {{
            const input = document.getElementById('searchInput');
            const filter = input.value.toUpperCase();
            const table = document.getElementById('samplesTable');
            const tr = table.getElementsByTagName('tr');
            
            for (let i = 1; i < tr.length; i++) {{
                const td = tr[i].getElementsByTagName('td')[0];
                if (td) {{
                    const txtValue = td.textContent || td.innerText;
                    tr[i].style.display = txtValue.toUpperCase().indexOf(filter) > -1 ? '' : 'none';
                }}
            }}
        }}
        
        function sortTable(columnIndex) {{
            const table = document.getElementById('samplesTable');
            const tbody = table.querySelector('tbody');
            const rows = Array.from(tbody.querySelectorAll('tr'));
            
            rows.sort((a, b) => {{
                const aText = a.cells[columnIndex].textContent.trim();
                const bText = b.cells[columnIndex].textContent.trim();
                return aText.localeCompare(bText, undefined, {{numeric: true}});
            }});
            
            rows.forEach(row => tbody.appendChild(row));
        }}
    </script>
</head>
<body>
    <div class="dashboard">
        <div class="header">
            <h1>🧬 SMN Batch Analysis Dashboard</h1>
            <p>Comprehensive SMA screening results</p>
            <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-number info">{total_samples}</div>
                <div class="stat-label">Total Samples</div>
            </div>
            <div class="stat-card">
                <div class="stat-number critical">{status_counts.get('AFFECTED', 0)}</div>
                <div class="stat-label">SMA Affected</div>
            </div>
            <div class="stat-card">
                <div class="stat-number warning">{status_counts.get('CARRIER', 0)}</div>
                <div class="stat-label">SMA Carriers</div>
            </div>
            <div class="stat-card">
                <div class="stat-number success">{status_counts.get('NORMAL', 0)}</div>
                <div class="stat-label">Normal Results</div>
            </div>
        </div>
        
        <div class="section">
            <h3>📊 Batch Statistics</h3>
            {stats_html}
        </div>
        
        <div class="section">
            <h3>🔍 Sample Details</h3>
            <input type="text" id="searchInput" onkeyup="filterTable()" 
                   placeholder="Search samples..." style="margin-bottom: 15px; padding: 8px; width: 300px;">
            {table_html}
        </div>
        
        {plots_html}
        
        <div class="footer">
            <p><strong>Pipeline:</strong> SMN CNV Detection v2.0-SMA-OPTIMIZED</p>
            <p><strong>Analysis Focus:</strong> SMN1/SMN2 Exons 7 & 8</p>
            <p><strong>Disclaimer:</strong> For research and clinical decision support. Confirmatory testing recommended.</p>
        </div>
    </div>
</body>
</html>
"""
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        logging.info(f"Batch HTML dashboard saved: {output_file}")
    
    def _generate_statistics_table(self, df):
        """Generate statistics table HTML"""
        stats = {
            'Total Samples': len(df),
            'Reference Samples': len(df[df['sample_type'] == 'reference']),
            'Test Samples': len(df[df['sample_type'] == 'test']),
            'Quality PASS': len(df[df['quality_flag'] == 'PASS']),
            'Quality WARNING': len(df[df['quality_flag'] == 'WARNING']),
            'Quality FAIL': len(df[df['quality_flag'] == 'FAIL']),
            'Average Quality Score': f"{df['quality_score'].mean():.3f}",
            'Samples Requiring Follow-up': len(df[
                df['sma_status'].isin(['AFFECTED', 'CARRIER', 'UNCERTAIN']) |
                df['quality_flag'].isin(['WARNING', 'FAIL'])
            ])
        }
        
        html = "<table><tr><th>Metric</th><th>Value</th></tr>"
        for metric, value in stats.items():
            html += f"<tr><td>{metric}</td><td>{value}</td></tr>"
        html += "</table>"
        
        return html
    
    def _generate_sample_table(self, df):
        """Generate sample details table HTML"""
        html = """
        <table id="samplesTable">
            <thead>
                <tr>
                    <th onclick="sortTable(0)" style="cursor: pointer;">Sample ID ↕</th>
                    <th onclick="sortTable(1)" style="cursor: pointer;">SMA Status ↕</th>
                    <th onclick="sortTable(2)" style="cursor: pointer;">Quality ↕</th>
                    <th onclick="sortTable(3)" style="cursor: pointer;">SMN1 CN ↕</th>
                    <th onclick="sortTable(4)" style="cursor: pointer;">SMN2 CN ↕</th>
                    <th>Actions</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for _, row in df.iterrows():
            sample_id = row['sample_id']
            sma_status = row.get('sma_status', 'UNKNOWN')
            quality_flag = row.get('quality_flag', 'UNKNOWN')
            smn1_cn = row.get('smn1_average_cn', 'N/A')
            smn2_cn = row.get('smn2_average_cn', 'N/A')
            
            html += f"""
                <tr>
                    <td>{sample_id}</td>
                    <td><span class="status-{sma_status.lower()}">{sma_status}</span></td>
                    <td><span class="quality-{quality_flag.lower()}">{quality_flag}</span></td>
                    <td>{smn1_cn}</td>
                    <td>{smn2_cn}</td>
                    <td><a href="{sample_id}/{sample_id}_report.html" target="_blank">View Report</a></td>
                </tr>
            """
        
        html += "</tbody></table>"
        return html
    
    def _generate_plots_html(self, df):
        """Generate plots section for HTML dashboard"""
        # This would generate matplotlib plots and embed them as base64
        # For simplicity, returning placeholder
        return """
        <div class="section">
            <h3>📈 Analysis Visualizations</h3>
            <div class="plots">
                <p>Interactive plots would be generated here showing:</p>
                <ul>
                    <li>SMA status distribution</li>
                    <li>Quality score distribution</li>
                    <li>Copy number scatter plots</li>
                    <li>Threshold performance metrics</li>
                </ul>
            </div>
        </div>
        """

def main():
    parser = argparse.ArgumentParser(description='Generate batch analysis summaries')
    parser.add_argument('--reports-dir', required=True, help='Directory containing sample reports')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--include-mlpa-concordance', action='store_true',
                       help='Include MLPA concordance analysis')
    parser.add_argument('--threshold-performance-analysis', action='store_true',
                       help='Include threshold performance analysis')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Initialize generator
    generator = BatchSummaryGenerator()
    
    # Load sample reports
    logger.info("Loading sample report data")
    df = generator.load_sample_reports(args.reports_dir)
    
    if df.empty:
        logger.error("No sample reports found")
        sys.exit(1)
    
    logger.info(f"Loaded data for {len(df)} samples")
    
    # Generate TSV summary
    tsv_file = Path(args.output_dir) / f"batch_summary_{generator.timestamp}.tsv"
    summary_df = generator.generate_batch_tsv(df, tsv_file)
    
    # Generate HTML dashboard
    html_file = Path(args.output_dir) / f"batch_summary_{generator.timestamp}.html"
    generator.generate_batch_html_dashboard(df, html_file)
    
    # Log key findings
    critical_samples = len(df[df['sma_status'] == 'AFFECTED'])
    carrier_samples = len(df[df['sma_status'] == 'CARRIER'])
    failed_samples = len(df[df['quality_flag'] == 'FAIL'])
    
    logger.info("Batch summary generation completed:")
    logger.info(f"  Critical findings (AFFECTED): {critical_samples}")
    logger.info(f"  Carrier samples: {carrier_samples}")
    logger.info(f"  Failed quality control: {failed_samples}")
    logger.info(f"  TSV summary: {tsv_file}")
    logger.info(f"  HTML dashboard: {html_file}")

if __name__ == '__main__':
    main()
