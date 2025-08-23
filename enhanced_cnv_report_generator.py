#!/usr/bin/env python3

"""
Enhanced SMN CNV Report Generator with Population Evidence
Generates comprehensive TXT and HTML reports with clinical context
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from dataclasses import dataclass
import logging

@dataclass
class SMNClinicalContext:
    """Clinical context from population databases"""
    clinvar_pathogenic_variants: int
    population_frequencies: Dict
    literature_links: Dict
    clinical_guidelines: Dict
    sma_phenotype_correlation: Dict

class SMNReportGenerator:
    """Generates comprehensive SMN CNV reports with population evidence"""
    
    def __init__(self, results_dir: str, population_cache: str):
        self.results_dir = Path(results_dir)
        self.population_cache = Path(population_cache)
        self.template_dir = self.results_dir / "templates"
        self.template_dir.mkdir(exist_ok=True)
        
        # Clinical interpretation mappings
        self.sma_risk_mapping = {
            (0, 0): {"risk": "AFFECTED", "description": "SMN1 homozygous deletion - Likely SMA affected"},
            (1, 0): {"risk": "AFFECTED", "description": "SMN1 deletion with SMN2 loss - Severe SMA risk"},
            (0, 1): {"risk": "AFFECTED", "description": "SMN1 loss with minimal SMN2 - High SMA risk"},
            (1, 1): {"risk": "CARRIER", "description": "SMN1 heterozygous deletion - SMA carrier"},
            (1, 2): {"risk": "CARRIER", "description": "SMN1 carrier with normal SMN2"},
            (2, 0): {"risk": "MODIFIED", "description": "Normal SMN1, SMN2 deletion - Atypical pattern"},
            (2, 1): {"risk": "LOW", "description": "Normal SMN1, reduced SMN2 - Low risk"},
            (2, 2): {"risk": "NORMAL", "description": "Normal SMN1 and SMN2 copy numbers"},
            (3, 2): {"risk": "LOW", "description": "SMN1 duplication may be protective"},
            (2, 3): {"risk": "LOW", "description": "SMN2 duplication may modify phenotype"}
        }
        
    def generate_sample_report(self, sample_result: Dict, output_dir: Path) -> Dict:
        """Generate comprehensive sample report"""
        sample_id = sample_result['sample_id']
        sample_dir = output_dir / sample_id
        sample_dir.mkdir(exist_ok=True)
        
        # Load population evidence
        population_context = self._load_population_context()
        
        # Perform clinical interpretation
        clinical_assessment = self._assess_clinical_risk(sample_result, population_context)
        
        # Generate reports
        txt_report = self._generate_txt_report(sample_result, clinical_assessment, population_context)
        html_report = self._generate_html_report(sample_result, clinical_assessment, population_context)
        
        # Save reports
        txt_file = sample_dir / f"{sample_id}_report.txt"
        html_file = sample_dir / f"{sample_id}_report.html"
        json_file = sample_dir / f"{sample_id}_data.json"
        
        with open(txt_file, 'w') as f:
            f.write(txt_report)
        
        with open(html_file, 'w') as f:
            f.write(html_report)
        
        # Save structured data
        report_data = {
            'sample_id': sample_id,
            'analysis_timestamp': datetime.now().isoformat(),
            'smn_results': sample_result,
            'clinical_assessment': clinical_assessment,
            'population_context': population_context,
            'report_version': '2.0'
        }
        
        with open(json_file, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        return {
            'sample_id': sample_id,
            'txt_report': txt_file,
            'html_report': html_file,
            'json_data': json_file,
            'clinical_risk': clinical_assessment['overall_risk'],
            'sma_category': clinical_assessment['sma_category']
        }
    
    def _load_population_context(self) -> Dict:
        """Load cached population evidence"""
        context = {}
        
        # Load ClinVar data
        clinvar_file = self.population_cache / "clinvar_SMN1.json"
        if clinvar_file.exists():
            with open(clinvar_file, 'r') as f:
                context['clinvar'] = json.load(f)
        
        # Load population frequencies
        pop_freq_file = self.population_cache / "population_frequencies.json"
        if pop_freq_file.exists():
            with open(pop_freq_file, 'r') as f:
                context['population_frequencies'] = json.load(f)
        
        # Static clinical guidelines and phenotype correlations
        context['clinical_guidelines'] = self._get_clinical_guidelines()
        context['phenotype_correlations'] = self._get_phenotype_correlations()
        context['literature_links'] = self._get_literature_links()
        
        return context
    
    def _get_clinical_guidelines(self) -> Dict:
        """Clinical guidelines for SMN CNV interpretation"""
        return {
            'acmg_recommendations': {
                'smn1_deletion_reporting': 'Report homozygous and heterozygous SMN1 deletions',
                'smn2_copy_modifier': 'SMN2 copy number modifies SMA severity',
                'carrier_counseling': 'Heterozygous SMN1 deletion carriers need genetic counseling'
            },
            'sma_severity_prediction': {
                'smn1_0_smn2_1-2': 'SMA Type I (severe)',
                'smn1_0_smn2_3-4': 'SMA Type II-III (moderate-mild)',
                'smn1_1_any_smn2': 'SMA carrier - reproductive counseling recommended'
            },
            'testing_recommendations': {
                'confirmatory_testing': 'MLPA or qPCR confirmation recommended',
                'family_screening': 'Cascade testing for at-risk family members',
                'reproductive_counseling': 'Preconception/prenatal counseling for carriers'
            }
        }
    
    def _get_phenotype_correlations(self) -> Dict:
        """SMN copy number to phenotype correlations"""
        return {
            'smn1_0_copies': {
                'phenotype': 'Spinal Muscular Atrophy (SMA)',
                'severity_modifiers': 'SMN2 copy number determines severity',
                'smn2_1_copy': 'SMA Type I (severe, onset <6 months)',
                'smn2_2_copies': 'SMA Type I-II (severe to intermediate)',
                'smn2_3_copies': 'SMA Type II-III (intermediate to mild)',
                'smn2_4_copies': 'SMA Type III-IV (mild to asymptomatic)'
            },
            'smn1_1_copy': {
                'phenotype': 'SMA carrier',
                'reproductive_risk': '25% chance of affected offspring if partner is carrier',
                'population_frequency': '1 in 40-60 individuals'
            },
            'smn1_2_copies': {
                'phenotype': 'Normal',
                'note': 'Standard copy number for most individuals'
            },
            'smn1_3_plus_copies': {
                'phenotype': 'Normal with potential protective effect',
                'note': 'May provide additional functional SMN protein'
            }
        }
    
    def _get_literature_links(self) -> Dict:
        """Generate literature and database links"""
        return {
            'primary_databases': {
                'clinvar_smn1': 'https://www.ncbi.nlm.nih.gov/clinvar/?term=SMN1%5Bgene%5D',
                'omim_sma': 'https://www.omim.org/entry/253300',
                'ncbi_smn1': 'https://www.ncbi.nlm.nih.gov/gene/6606',
                'ncbi_smn2': 'https://www.ncbi.nlm.nih.gov/gene/6607'
            },
            'literature_searches': {
                'pubmed_sma_cnv': 'https://pubmed.ncbi.nlm.nih.gov/?term=SMN1+copy+number+spinal+muscular+atrophy',
                'pubmed_carrier_screening': 'https://pubmed.ncbi.nlm.nih.gov/?term=SMN1+carrier+screening',
                'pubmed_sma_severity': 'https://pubmed.ncbi.nlm.nih.gov/?term=SMN2+copy+number+SMA+severity'
            },
            'clinical_resources': {
                'acmg_guidelines': 'https://www.acmg.net/PDFLibrary/SMA-Technical-Standard.pdf',
                'sma_foundation': 'https://www.curesma.org/',
                'genetests': 'https://www.ncbi.nlm.nih.gov/books/NBK1352/'
            }
        }
    
    def _assess_clinical_risk(self, sample_result: Dict, population_context: Dict) -> Dict:
        """Comprehensive clinical risk assessment"""
        # Extract copy numbers
        smn1_cn = self._get_gene_copy_number(sample_result, 'SMN1')
        smn2_cn = self._get_gene_copy_number(sample_result, 'SMN2')
        
        # Get clinical interpretation
        risk_key = (smn1_cn, smn2_cn)
        clinical_interpretation = self.sma_risk_mapping.get(
            risk_key, 
            {"risk": "UNCERTAIN", "description": f"Atypical SMN1={smn1_cn}, SMN2={smn2_cn} pattern"}
        )
        
        # Population frequency context
        pop_freq = population_context.get('population_frequencies', {})
        carrier_freq = pop_freq.get('smn1_carrier_frequency', {}).get('global', 'Unknown')
        
        # Quality control assessment
        qc_status = self._assess_quality_control(sample_result)
        
        # Evidence strength
        evidence_strength = self._assess_evidence_strength(sample_result)
        
        return {
            'smn1_copy_number': smn1_cn,
            'smn2_copy_number': smn2_cn,
            'overall_risk': clinical_interpretation['risk'],
            'sma_category': clinical_interpretation['description'],
            'population_frequency': carrier_freq,
            'quality_control': qc_status,
            'evidence_strength': evidence_strength,
            'recommendations': self._generate_recommendations(clinical_interpretation['risk'], qc_status),
            'interpretation_confidence': self._calculate_confidence(sample_result, evidence_strength)
        }
    
    def _get_gene_copy_number(self, sample_result: Dict, gene: str) -> int:
        """Extract gene-level copy number from sample results"""
        exon7_key = f"{gene}_exon7"
        exon8_key = f"{gene}_exon8"
        
        exon_depths = sample_result.get('exon_depths', {})
        
        # Get copy numbers from both exons
        exon7_cn = exon_depths.get(exon7_key, {}).get('copy_number', 2)
        exon8_cn = exon_depths.get(exon8_key, {}).get('copy_number', 2)
        
        # Use median of available exons
        valid_cns = [cn for cn in [exon7_cn, exon8_cn] if cn is not None and not np.isnan(cn)]
        
        if valid_cns:
            return int(np.median(valid_cns))
        else:
            return 2  # Default to normal if no data
    
    def _assess_quality_control(self, sample_result: Dict) -> Dict:
        """Assess quality control metrics"""
        exon_depths = sample_result.get('exon_depths', {})
        fallback_used = sample_result.get('fallback_used', False)
        
        qc_flags = []
        overall_quality = "PASS"
        
        # Check for adequate read coverage
        low_coverage_exons = []
        for exon_id, depth_data in exon_depths.items():
            if not depth_data.get('reads_detected', True):
                qc_flags.append(f"No reads detected in {exon_id}")
                low_coverage_exons.append(exon_id)
            elif depth_data.get('mean_depth', 0) < 20:
                qc_flags.append(f"Low coverage in {exon_id} ({depth_data.get('mean_depth', 0):.1f}x)")
                low_coverage_exons.append(exon_id)
        
        if low_coverage_exons:
            if len(low_coverage_exons) >= 3:
                overall_quality = "FAIL"
            else:
                overall_quality = "WARNING"
        
        if fallback_used:
            qc_flags.append("Fallback analysis used due to insufficient primary exon coverage")
            overall_quality = "WARNING" if overall_quality == "PASS" else overall_quality
        
        return {
            'overall_quality': overall_quality,
            'qc_flags': qc_flags,
            'low_coverage_exons': low_coverage_exons,
            'fallback_used': fallback_used,
            'fallback_regions': sample_result.get('fallback_regions', [])
        }
    
    def _assess_evidence_strength(self, sample_result: Dict) -> str:
        """Assess strength of evidence for CNV call"""
        exon_depths = sample_result.get('exon_depths', {})
        
        # Count exons with strong evidence
        strong_evidence_count = 0
        total_exons = 0
        
        for exon_id, depth_data in exon_depths.items():
            if exon_id.endswith(('exon7', 'exon8')):  # Focus on critical exons
                total_exons += 1
                
                if (depth_data.get('reads_detected', False) and 
                    depth_data.get('mean_depth', 0) >= 20 and
                    depth_data.get('confidence', 'low') in ['high', 'medium']):
                    strong_evidence_count += 1
        
        if strong_evidence_count == total_exons and total_exons >= 4:
            return "STRONG"
        elif strong_evidence_count >= total_exons * 0.75:
            return "MODERATE"
        elif strong_evidence_count >= total_exons * 0.5:
            return "LIMITED"
        else:
            return "WEAK"
    
    def _calculate_confidence(self, sample_result: Dict, evidence_strength: str) -> str:
        """Calculate overall confidence in interpretation"""
        qc_status = self._assess_quality_control(sample_result)['overall_quality']
        
        if qc_status == "PASS" and evidence_strength == "STRONG":
            return "HIGH"
        elif qc_status in ["PASS", "WARNING"] and evidence_strength in ["STRONG", "MODERATE"]:
            return "MEDIUM"
        elif qc_status != "FAIL" and evidence_strength in ["MODERATE", "LIMITED"]:
            return "LOW"
        else:
            return "VERY_LOW"
    
    def _generate_recommendations(self, risk_category: str, qc_status: Dict) -> List[str]:
        """Generate clinical recommendations"""
        recommendations = []
        
        # Risk-based recommendations
        if risk_category == "AFFECTED":
            recommendations.extend([
                "Confirmatory testing with MLPA or qPCR is strongly recommended",
                "Refer to medical genetics for comprehensive evaluation",
                "Consider neurological evaluation for SMA signs/symptoms",
                "Genetic counseling for family planning"
            ])
        elif risk_category == "CARRIER":
            recommendations.extend([
                "Genetic counseling recommended",
                "Partner screening advised for family planning",
                "Confirmatory testing may be considered",
                "Cascade testing for at-risk family members"
            ])
        elif risk_category == "NORMAL":
            recommendations.append("No immediate clinical action required based on CNV analysis")
        elif risk_category in ["UNCERTAIN", "MODIFIED"]:
            recommendations.extend([
                "Confirmatory testing recommended due to atypical pattern",
                "Clinical correlation and genetic counseling advised",
                "Consider additional molecular testing"
            ])
        
        # Quality-based recommendations
        if qc_status['overall_quality'] == "FAIL":
            recommendations.insert(0, "CRITICAL: Poor sample quality - repeat testing strongly recommended")
        elif qc_status['overall_quality'] == "WARNING":
            recommendations.append("Consider repeat testing due to quality concerns")
        
        if qc_status['fallback_used']:
            recommendations.append("Confirmatory testing recommended due to limited primary exon coverage")
        
        return recommendations
    
    def _generate_txt_report(self, sample_result: Dict, clinical_assessment: Dict, population_context: Dict) -> str:
        """Generate detailed TXT report"""
        sample_id = sample_result['sample_id']
        
        report = f"""
SMN CNV ANALYSIS REPORT
=======================

Sample ID: {sample_id}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Report Version: 2.0 (Optimized for SMA Detection)

CLINICAL SUMMARY
================
SMN1 Copy Number: {clinical_assessment['smn1_copy_number']}
SMN2 Copy Number: {clinical_assessment['smn2_copy_number']}
SMA Risk Category: {clinical_assessment['overall_risk']}
Clinical Interpretation: {clinical_assessment['sma_category']}
Interpretation Confidence: {clinical_assessment['interpretation_confidence']}

QUALITY CONTROL
===============
Overall Quality: {clinical_assessment['quality_control']['overall_quality']}
Evidence Strength: {clinical_assessment['evidence_strength']}
"""
        
        # Add QC flags if present
        qc_flags = clinical_assessment['quality_control']['qc_flags']
        if qc_flags:
            report += "\nQuality Control Flags:\n"
            for flag in qc_flags:
                report += f"- {flag}\n"
        
        # Add fallback information
        if clinical_assessment['quality_control']['fallback_used']:
            report += f"\nFallback Analysis Used:\n"
            for region in clinical_assessment['quality_control']['fallback_regions']:
                report += f"- {region}\n"
        
        # Detailed exon results
        report += "\nDETAILED EXON RESULTS\n"
        report += "=====================\n"
        
        exon_depths = sample_result.get('exon_depths', {})
        for exon_id in ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']:
            if exon_id in exon_depths:
                depth_data = exon_depths[exon_id]
                report += f"\n{exon_id}:\n"
                report += f"  Reads Detected: {'Yes' if depth_data.get('reads_detected', False) else 'NO - FALLBACK ANALYSIS'}\n"
                report += f"  Mean Depth: {depth_data.get('mean_depth', 0):.1f}x\n"
                report += f"  Z-score: {depth_data.get('z_score', 'N/A')}\n"
                report += f"  Copy Number: {depth_data.get('copy_number', 'N/A')}\n"
                report += f"  Confidence: {depth_data.get('confidence', 'N/A')}\n"
        
        # Clinical recommendations
        report += "\nCLINICAL RECOMMENDATIONS\n"
        report += "========================\n"
        for i, rec in enumerate(clinical_assessment['recommendations'], 1):
            report += f"{i}. {rec}\n"
        
        # Population context
        pop_freq = population_context.get('population_frequencies', {})
        if pop_freq:
            report += "\nPOPULATION CONTEXT\n"
            report += "==================\n"
            
            carrier_freq = pop_freq.get('smn1_carrier_frequency', {})
            if carrier_freq:
                report += "SMN1 Carrier Frequencies:\n"
                for population, frequency in carrier_freq.items():
                    if population != 'global':
                        report += f"  {population.title()}: {frequency}\n"
            
            sma_incidence = pop_freq.get('sma_incidence', {})
            if sma_incidence:
                report += f"\nSMA Disease Incidence: {sma_incidence.get('global', 'Unknown')}\n"
        
        # Evidence sources
        report += "\nEVIDENCE SOURCES & REFERENCES\n"
        report += "=============================\n"
        
        clinvar_data = population_context.get('clinvar', {})
        if clinvar_data:
            report += f"ClinVar Pathogenic SMN1 Variants: {clinvar_data.get('pathogenic_variants', 'N/A')}\n"
        
        literature_links = population_context.get('literature_links', {})
        if literature_links:
            report += "\nRelevant Resources:\n"
            for category, links in literature_links.items():
                report += f"\n{category.replace('_', ' ').title()}:\n"
                for name, url in links.items():
                    report += f"  {name.replace('_', ' ').title()}: {url}\n"
        
        # Technical details
        report += "\nTECHNICAL DETAILS\n"
        report += "=================\n"
        report += "Analysis Focus: SMN1/SMN2 exons 7 and 8 (critical SMA loci)\n"
        report += "Normalization: Z-score based with MLPA-trained thresholds\n"
        report += "Copy Number Calling: Adaptive thresholds with population validation\n"
        report += "Quality Metrics: Read depth, coverage uniformity, Z-score confidence\n"
        
        report += "\nDISCLAIMER\n"
        report += "==========\n"
        report += "This analysis is for research purposes. Clinical decisions should not be made\n"
        report += "solely based on this report. Confirmatory testing and genetic counseling\n"
        report += "are recommended for actionable findings.\n"
        
        return report
    
    def _generate_html_report(self, sample_result: Dict, clinical_assessment: Dict, population_context: Dict) -> str:
        """Generate comprehensive HTML report"""
        sample_id = sample_result['sample_id']
        risk_category = clinical_assessment['overall_risk']
        
        # Determine alert class based on risk
        alert_class = {
            'AFFECTED': 'alert-danger',
            'CARRIER': 'alert-warning', 
            'NORMAL': 'alert-success',
            'UNCERTAIN': 'alert-info',
            'MODIFIED': 'alert-info'
        }.get(risk_category, 'alert-secondary')
        
        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN CNV Analysis Report - {sample_id}</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        .risk-badge {{
            font-size: 1.2em;
            padding: 0.5em 1em;
        }}
        .exon-result {{
            border-left: 4px solid #007bff;
            padding: 1rem;
            margin: 0.5rem 0;
            background-color: #f8f9fa;
        }}
        .no-reads {{
            border-left-color: #dc3545;
            background-color: #f8d7da;
        }}
        .evidence-section {{
            background-color: #e3f2fd;
            border: 1px solid #2196f3;
            border-radius: 0.5rem;
            padding: 1rem;
            margin: 1rem 0;
        }}
        .population-data {{
            background-color: #f3e5f5;
            border: 1px solid #9c27b0;
            border-radius: 0.5rem;
            padding: 1rem;
            margin: 1rem 0;
        }}
        .footer {{
            background-color: #f8f9fa;
            padding: 2rem;
            margin-top: 3rem;
            border-top: 1px solid #dee2e6;
        }}
    </style>
</head>
<body>
    <div class="container-fluid py-4">
        <div class="row">
            <div class="col-lg-12">
                
                <!-- Header -->
                <div class="card mb-4">
                    <div class="card-header bg-primary text-white">
                        <h1 class="mb-0"><i class="fas fa-dna me-2"></i>SMN CNV Analysis Report</h1>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h3>Sample: <strong>{sample_id}</strong></h3>
                                <p class="text-muted">Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                                <p class="text-muted">Report Version: 2.0 (SMA-Optimized)</p>
                            </div>
                            <div class="col-md-6 text-end">
                                <div class="alert {alert_class} risk-badge">
                                    <strong>SMA Risk: {risk_category}</strong>
                                </div>
                                <p class="text-muted">Confidence: {clinical_assessment['interpretation_confidence']}</p>
                            </div>
                        </div>
                    </div>
                </div>
                
                <!-- Clinical Summary -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-stethoscope me-2"></i>Clinical Summary</h2>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h4>Copy Number Results</h4>
                                <table class="table table-bordered">
                                    <tr>
                                        <td><strong>SMN1 Copy Number</strong></td>
                                        <td class="text-center">
                                            <span class="badge bg-primary fs-6">{clinical_assessment['smn1_copy_number']}</span>
                                        </td>
                                    </tr>
                                    <tr>
                                        <td><strong>SMN2 Copy Number</strong></td>
                                        <td class="text-center">
                                            <span class="badge bg-info fs-6">{clinical_assessment['smn2_copy_number']}</span>
                                        </td>
                                    </tr>
                                </table>
                            </div>
                            <div class="col-md-6">
                                <h4>Clinical Interpretation</h4>
                                <div class="alert {alert_class}">
                                    <strong>{clinical_assessment['sma_category']}</strong>
                                </div>
                                <p><strong>Evidence Strength:</strong> {clinical_assessment['evidence_strength']}</p>
                                <p><strong>Population Context:</strong> {clinical_assessment['population_frequency']} carrier frequency</p>
                            </div>
                        </div>
                    </div>
                </div>
        """
        
        # Quality Control Section
        qc_status = clinical_assessment['quality_control']
        qc_alert_class = {
            'PASS': 'alert-success',
            'WARNING': 'alert-warning',
            'FAIL': 'alert-danger'
        }.get(qc_status['overall_quality'], 'alert-secondary')
        
        html += f"""
                <!-- Quality Control -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-check-circle me-2"></i>Quality Control</h2>
                    </div>
                    <div class="card-body">
                        <div class="alert {qc_alert_class}">
                            <strong>Overall Quality: {qc_status['overall_quality']}</strong>
                        </div>
        """
        
        if qc_status['qc_flags']:
            html += "<h5>Quality Control Flags:</h5><ul>"
            for flag in qc_status['qc_flags']:
                html += f"<li>{flag}</li>"
            html += "</ul>"
        
        if qc_status['fallback_used']:
            html += f"""
                        <div class="alert alert-info">
                            <strong>Fallback Analysis Used</strong><br>
                            Primary exon coverage was insufficient. Analysis used:
                            <ul>
            """
            for region in qc_status['fallback_regions']:
                html += f"<li>{region}</li>"
            html += "</ul></div>"
        
        html += """
                    </div>
                </div>
        """
        
        # Detailed Exon Results
        html += """
                <!-- Detailed Results -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-chart-bar me-2"></i>Detailed Exon Results</h2>
                    </div>
                    <div class="card-body">
        """
        
        exon_depths = sample_result.get('exon_depths', {})
        for exon_id in ['SMN1_exon7', 'SMN1_exon8', 'SMN2_exon7', 'SMN2_exon8']:
            if exon_id in exon_depths:
                depth_data = exon_depths[exon_id]
                reads_detected = depth_data.get('reads_detected', False)
                exon_class = "exon-result" + (" no-reads" if not reads_detected else "")
                
                html += f"""
                        <div class="{exon_class}">
                            <div class="row">
                                <div class="col-md-3">
                                    <h5>{exon_id.replace('_', ' ')}</h5>
                                    {'<i class="fas fa-exclamation-triangle text-danger"></i> No Reads' if not reads_detected else '<i class="fas fa-check-circle text-success"></i> Reads Detected'}
                                </div>
                                <div class="col-md-9">
                                    <div class="row">
                                        <div class="col-sm-3">
                                            <strong>Depth:</strong> {depth_data.get('mean_depth', 0):.1f}x
                                        </div>
                                        <div class="col-sm-3">
                                            <strong>Z-score:</strong> {depth_data.get('z_score', 'N/A')}
                                        </div>
                                        <div class="col-sm-3">
                                            <strong>Copy Number:</strong> {depth_data.get('copy_number', 'N/A')}
                                        </div>
                                        <div class="col-sm-3">
                                            <strong>Confidence:</strong> {depth_data.get('confidence', 'N/A')}
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                """
        
        html += """
                    </div>
                </div>
        """
        
        # Clinical Recommendations
        html += f"""
                <!-- Recommendations -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-clipboard-list me-2"></i>Clinical Recommendations</h2>
                    </div>
                    <div class="card-body">
                        <ol>
        """
        
        for rec in clinical_assessment['recommendations']:
            html += f"<li>{rec}</li>"
        
        html += """
                        </ol>
                    </div>
                </div>
        """
        
        # Population Evidence
        pop_freq = population_context.get('population_frequencies', {})
        clinvar_data = population_context.get('clinvar', {})
        
        html += """
                <!-- Population Evidence -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-globe me-2"></i>Population Evidence & Context</h2>
                    </div>
                    <div class="card-body">
        """
        
        if pop_freq:
            html += """
                        <div class="population-data">
                            <h4><i class="fas fa-users me-2"></i>Population Frequencies</h4>
            """
            
            carrier_freq = pop_freq.get('smn1_carrier_frequency', {})
            if carrier_freq:
                html += "<h5>SMN1 Carrier Frequencies by Population:</h5><ul>"
                for population, frequency in carrier_freq.items():
                    html += f"<li><strong>{population.title()}:</strong> {frequency}</li>"
                html += "</ul>"
            
            sma_incidence = pop_freq.get('sma_incidence', {})
            if sma_incidence:
                html += f"<p><strong>SMA Disease Incidence:</strong> {sma_incidence.get('global', 'Unknown')}</p>"
            
            html += "</div>"
        
        if clinvar_data:
            html += f"""
                        <div class="evidence-section">
                            <h4><i class="fas fa-database me-2"></i>ClinVar Evidence</h4>
                            <p><strong>Pathogenic SMN1 Variants:</strong> {clinvar_data.get('pathogenic_variants', 'N/A')}</p>
                            <p><a href="{clinvar_data.get('query_url', '#')}" target="_blank" class="btn btn-outline-primary btn-sm">
                                <i class="fas fa-external-link-alt"></i> View in ClinVar
                            </a></p>
                        </div>
            """
        
        # Literature Links
        literature_links = population_context.get('literature_links', {})
        if literature_links:
            html += """
                        <div class="evidence-section">
                            <h4><i class="fas fa-book-open me-2"></i>Literature & Resources</h4>
                            <div class="row">
            """
            
            for category, links in literature_links.items():
                html += f"""
                                <div class="col-md-4 mb-3">
                                    <h6>{category.replace('_', ' ').title()}</h6>
                                    <ul class="list-unstyled">
                """
                for name, url in links.items():
                    html += f"""
                                        <li><a href="{url}" target="_blank" class="text-decoration-none">
                                            <i class="fas fa-external-link-alt me-1"></i>{name.replace('_', ' ').title()}
                                        </a></li>
                    """
                html += "</ul></div>"
            
            html += """
                            </div>
                        </div>
            """
        
        html += """
                    </div>
                </div>
        """
        
        # Phenotype Correlations
        phenotype_data = population_context.get('phenotype_correlations', {})
        if phenotype_data:
            html += """
                <!-- Phenotype Correlations -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-user-md me-2"></i>Phenotype Correlations</h2>
                    </div>
                    <div class="card-body">
            """
            
            # Show relevant phenotype information based on copy numbers
            smn1_cn = clinical_assessment['smn1_copy_number']
            
            if smn1_cn == 0 and 'smn1_0_copies' in phenotype_data:
                smn1_0_data = phenotype_data['smn1_0_copies']
                html += f"""
                        <div class="alert alert-warning">
                            <h5>SMN1 Homozygous Deletion (0 copies)</h5>
                            <p><strong>Phenotype:</strong> {smn1_0_data.get('phenotype', 'N/A')}</p>
                            <p><strong>Severity Modifier:</strong> {smn1_0_data.get('severity_modifiers', 'N/A')}</p>
                            <h6>SMN2 Copy Number Effects:</h6>
                            <ul>
                                <li><strong>1 SMN2 copy:</strong> {smn1_0_data.get('smn2_1_copy', 'N/A')}</li>
                                <li><strong>2 SMN2 copies:</strong> {smn1_0_data.get('smn2_2_copies', 'N/A')}</li>
                                <li><strong>3 SMN2 copies:</strong> {smn1_0_data.get('smn2_3_copies', 'N/A')}</li>
                                <li><strong>4+ SMN2 copies:</strong> {smn1_0_data.get('smn2_4_copies', 'N/A')}</li>
                            </ul>
                        </div>
                """
            elif smn1_cn == 1 and 'smn1_1_copy' in phenotype_data:
                smn1_1_data = phenotype_data['smn1_1_copy']
                html += f"""
                        <div class="alert alert-info">
                            <h5>SMN1 Heterozygous Deletion (1 copy)</h5>
                            <p><strong>Phenotype:</strong> {smn1_1_data.get('phenotype', 'N/A')}</p>
                            <p><strong>Reproductive Risk:</strong> {smn1_1_data.get('reproductive_risk', 'N/A')}</p>
                            <p><strong>Population Frequency:</strong> {smn1_1_data.get('population_frequency', 'N/A')}</p>
                        </div>
                """
            elif smn1_cn == 2 and 'smn1_2_copies' in phenotype_data:
                smn1_2_data = phenotype_data['smn1_2_copies']
                html += f"""
                        <div class="alert alert-success">
                            <h5>SMN1 Normal (2 copies)</h5>
                            <p><strong>Phenotype:</strong> {smn1_2_data.get('phenotype', 'N/A')}</p>
                            <p><strong>Note:</strong> {smn1_2_data.get('note', 'N/A')}</p>
                        </div>
                """
            elif smn1_cn >= 3 and 'smn1_3_plus_copies' in phenotype_data:
                smn1_3_data = phenotype_data['smn1_3_plus_copies']
                html += f"""
                        <div class="alert alert-success">
                            <h5>SMN1 Duplication (3+ copies)</h5>
                            <p><strong>Phenotype:</strong> {smn1_3_data.get('phenotype', 'N/A')}</p>
                            <p><strong>Note:</strong> {smn1_3_data.get('note', 'N/A')}</p>
                        </div>
                """
            
            html += """
                    </div>
                </div>
            """
        
        # Technical Details
        html += """
                <!-- Technical Details -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h2 class="mb-0"><i class="fas fa-cogs me-2"></i>Technical Details</h2>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h5>Analysis Parameters</h5>
                                <ul>
                                    <li><strong>Target Regions:</strong> SMN1/SMN2 exons 7 and 8</li>
                                    <li><strong>Normalization:</strong> Z-score based with MLPA-trained thresholds</li>
                                    <li><strong>Copy Number Calling:</strong> Adaptive thresholds with population validation</li>
                                    <li><strong>Quality Metrics:</strong> Read depth, coverage uniformity, Z-score confidence</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h5>Clinical Focus</h5>
                                <ul>
                                    <li><strong>Primary Indication:</strong> Spinal Muscular Atrophy (SMA) screening</li>
                                    <li><strong>Critical Loci:</strong> SMN1 exons 7 & 8 (essential for SMA diagnosis)</li>
                                    <li><strong>Modifier Analysis:</strong> SMN2 copy number (severity prediction)</li>
                                    <li><strong>Population Context:</strong> Ancestry-specific carrier frequencies</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </div>
        """
        
        # Footer with disclaimer
        html += f"""
                <!-- Footer -->
                <div class="footer">
                    <div class="row">
                        <div class="col-md-8">
                            <h5><i class="fas fa-exclamation-triangle me-2 text-warning"></i>Important Disclaimer</h5>
                            <p class="mb-2">
                                This analysis is for research and clinical decision support purposes. 
                                Clinical decisions should not be made solely based on this computational analysis. 
                                Confirmatory testing with validated clinical methods (MLPA, qPCR) and genetic 
                                counseling are strongly recommended for actionable findings.
                            </p>
                            <p class="mb-0 text-muted">
                                <small>
                                    Report generated by SMN CNV Detection Pipeline v2.0 on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}. 
                                    This report integrates population evidence from ClinVar, population frequency databases, 
                                    and peer-reviewed literature to provide comprehensive clinical context.
                                </small>
                            </p>
                        </div>
                        <div class="col-md-4 text-end">
                            <div class="alert alert-light">
                                <strong>For Clinical Support:</strong><br>
                                <small>
                                    • Medical Genetics Consultation<br>
                                    • Genetic Counseling Services<br>
                                    • Confirmatory Laboratory Testing<br>
                                    • SMA Foundation Resources
                                </small>
                            </div>
                        </div>
                    </div>
                </div>
                
            </div>
        </div>
    </div>
    
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    <script>
        // Add interactive features
        document.addEventListener('DOMContentLoaded', function() {{
            // Add tooltips for technical terms
            var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
            var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {{
                return new bootstrap.Tooltip(tooltipTriggerEl);
            }});
            
            // Highlight critical findings
            const riskCategory = '{risk_category}';
            if (riskCategory === 'AFFECTED') {{
                document.body.style.borderTop = '5px solid #dc3545';
            }} else if (riskCategory === 'CARRIER') {{
                document.body.style.borderTop = '5px solid #ffc107';
            }}
        }});
    </script>
</body>
</html>
        """
        
        return html
    
    def generate_batch_report(self, all_sample_results: List[Dict], output_file: Path) -> Dict:
        """Generate consolidated batch report across all samples"""
        
        # Create batch summary data
        batch_data = {
            'total_samples': len(all_sample_results),
            'analysis_timestamp': datetime.now().isoformat(),
            'risk_distribution': {},
            'quality_distribution': {},
            'mlpa_concordance': None,
            'threshold_performance': {}
        }
        
        # Analyze risk distribution
        risk_counts = {}
        quality_counts = {}
        
        for result in all_sample_results:
            risk = result.get('clinical_risk', 'UNKNOWN')
            quality = result.get('quality_status', 'UNKNOWN')
            
            risk_counts[risk] = risk_counts.get(risk, 0) + 1
            quality_counts[quality] = quality_counts.get(quality, 0) + 1
        
        batch_data['risk_distribution'] = risk_counts
        batch_data['quality_distribution'] = quality_counts
        
        # Generate batch TSV report
        batch_tsv = self._generate_batch_tsv(all_sample_results, batch_data)
        
        # Generate batch HTML summary
        batch_html = self._generate_batch_html(all_sample_results, batch_data)
        
        # Save reports
        tsv_file = output_file.parent / f"batch_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"
        html_file = output_file.parent / f"batch_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
        
        with open(tsv_file, 'w') as f:
            f.write(batch_tsv)
        
        with open(html_file, 'w') as f:
            f.write(batch_html)
        
        return {
            'batch_tsv': tsv_file,
            'batch_html': html_file,
            'summary_stats': batch_data
        }
    
    def _generate_batch_tsv(self, all_results: List[Dict], batch_data: Dict) -> str:
        """Generate batch TSV summary"""
        
        tsv_lines = [
            "# SMN CNV Analysis Batch Summary",
            f"# Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"# Total Samples: {batch_data['total_samples']}",
            "#",
            "sample_id\tsmn1_copy_number\tsmn2_copy_number\tsma_risk_category\tquality_status\tevidence_strength\tinterpretation_confidence\tfallback_used\trecommendations"
        ]
        
        for result in all_results:
            sample_data = result.get('smn_results', {})
            clinical = result.get('clinical_assessment', {})
            
            line = f"{result.get('sample_id', 'Unknown')}\t"
            line += f"{clinical.get('smn1_copy_number', 'N/A')}\t"
            line += f"{clinical.get('smn2_copy_number', 'N/A')}\t"
            line += f"{clinical.get('overall_risk', 'N/A')}\t"
            line += f"{clinical.get('quality_control', {}).get('overall_quality', 'N/A')}\t"
            line += f"{clinical.get('evidence_strength', 'N/A')}\t"
            line += f"{clinical.get('interpretation_confidence', 'N/A')}\t"
            line += f"{clinical.get('quality_control', {}).get('fallback_used', False)}\t"
            line += f"{'; '.join(clinical.get('recommendations', []))}"
            
            tsv_lines.append(line)
        
        return "\n".join(tsv_lines)
    
    def _generate_batch_html(self, all_results: List[Dict], batch_data: Dict) -> str:
        """Generate batch HTML summary"""
        
        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SMN CNV Batch Analysis Summary</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
</head>
<body>
    <div class="container-fluid py-4">
        <div class="row">
            <div class="col-12">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <h1><i class="fas fa-chart-bar me-2"></i>SMN CNV Batch Analysis Summary</h1>
                        <p class="mb-0">Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Total Samples: {batch_data['total_samples']}</p>
                    </div>
                    <div class="card-body">
                        
                        <!-- Summary Statistics -->
                        <div class="row mb-4">
                            <div class="col-md-6">
                                <h3>SMA Risk Distribution</h3>
                                <canvas id="riskChart" width="400" height="200"></canvas>
                            </div>
                            <div class="col-md-6">
                                <h3>Quality Control Distribution</h3>
                                <canvas id="qualityChart" width="400" height="200"></canvas>
                            </div>
                        </div>
                        
                        <!-- Detailed Results Table -->
                        <h3>Detailed Results</h3>
                        <div class="table-responsive">
                            <table class="table table-striped table-hover">
                                <thead class="table-dark">
                                    <tr>
                                        <th>Sample ID</th>
                                        <th>SMN1 CN</th>
                                        <th>SMN2 CN</th>
                                        <th>SMA Risk</th>
                                        <th>Quality</th>
                                        <th>Evidence</th>
                                        <th>Confidence</th>
                                        <th>Fallback</th>
                                        <th>Actions</th>
                                    </tr>
                                </thead>
                                <tbody>
        """
        
        for result in all_results:
            sample_id = result.get('sample_id', 'Unknown')
            clinical = result.get('clinical_assessment', {})
            risk = clinical.get('overall_risk', 'N/A')
            
            # Determine row class based on risk
            row_class = {
                'AFFECTED': 'table-danger',
                'CARRIER': 'table-warning',
                'NORMAL': 'table-success',
                'UNCERTAIN': 'table-info'
            }.get(risk, '')
            
            html += f"""
                                    <tr class="{row_class}">
                                        <td><strong>{sample_id}</strong></td>
                                        <td>{clinical.get('smn1_copy_number', 'N/A')}</td>
                                        <td>{clinical.get('smn2_copy_number', 'N/A')}</td>
                                        <td><span class="badge bg-{'danger' if risk=='AFFECTED' else 'warning' if risk=='CARRIER' else 'success'}">{risk}</span></td>
                                        <td>{clinical.get('quality_control', {}).get('overall_quality', 'N/A')}</td>
                                        <td>{clinical.get('evidence_strength', 'N/A')}</td>
                                        <td>{clinical.get('interpretation_confidence', 'N/A')}</td>
                                        <td>{'Yes' if clinical.get('quality_control', {}).get('fallback_used', False) else 'No'}</td>
                                        <td><a href="{sample_id}/{sample_id}_report.html" class="btn btn-sm btn-outline-primary">View Report</a></td>
                                    </tr>
            """
        
        html += """
                                </tbody>
                            </table>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    
    <script>
        // Risk Distribution Chart
        const riskData = """ + json.dumps(batch_data['risk_distribution']) + """;
        const riskChart = new Chart(document.getElementById('riskChart'), {
            type: 'doughnut',
            data: {
                labels: Object.keys(riskData),
                datasets: [{
                    data: Object.values(riskData),
                    backgroundColor: ['#dc3545', '#ffc107', '#28a745', '#17a2b8', '#6c757d']
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false
            }
        });
        
        // Quality Distribution Chart
        const qualityData = """ + json.dumps(batch_data['quality_distribution']) + """;
        const qualityChart = new Chart(document.getElementById('qualityChart'), {
            type: 'bar',
            data: {
                labels: Object.keys(qualityData),
                datasets: [{
                    data: Object.values(qualityData),
                    backgroundColor: ['#28a745', '#ffc107', '#dc3545']
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        display: false
                    }
                }
            }
        });
    </script>
</body>
</html>
        """
        
        return html

def main():
    """Example usage of the report generator"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate SMN CNV reports')
    parser.add_argument('--results-dir', required=True, help='Results directory')
    parser.add_argument('--sample-data', required=True, help='Sample analysis data file (JSON)')
    parser.add_argument('--output-dir', required=True, help='Report output directory')
    
    args = parser.parse_args()
    
    # Initialize report generator
    generator = SMNReportGenerator(args.results_dir, f"{args.results_dir}/population_cache")
    
    # Load sample data
    with open(args.sample_data, 'r') as f:
        sample_data = json.load(f)
    
    # Generate report
    report_files = generator.generate_sample_report(sample_data, Path(args.output_dir))
    
    print(f"Generated reports:")
    print(f"  TXT: {report_files['txt_report']}")
    print(f"  HTML: {report_files['html_report']}")
    print(f"  JSON: {report_files['json_data']}")

if __name__ == "__main__":
    main()
