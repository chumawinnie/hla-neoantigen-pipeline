/*
========================================================================================
    GENERATE REPORT
========================================================================================
    Generate clinical summary report for molecular tumor board
*/

process GENERATE_REPORT {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/${meta.id}/report", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(filtered_neoantigens), path(ranked_neoantigens), path(neoantigen_summary), path(hla_alleles)

    output:
    tuple val(meta), path("${meta.id}_neoantigen_report.html"),     emit: html_report
    tuple val(meta), path("${meta.id}_neoantigen_report.pdf"),      emit: pdf_report, optional: true
    path "*.log",                                                    emit: logs

    script:
    """
    #!/usr/bin/env python3
    
    import json
    import pandas as pd
    from datetime import datetime
    from pathlib import Path
    
    # Load data
    with open('${neoantigen_summary}', 'r') as f:
        summary = json.load(f)
    
    hla_alleles = Path('${hla_alleles}').read_text().strip()
    
    # Load ranked neoantigens
    try:
        ranked_df = pd.read_csv('${ranked_neoantigens}', sep='\\t')
    except:
        ranked_df = pd.DataFrame()
    
    # Generate HTML report
    html_content = f'''<!DOCTYPE html>
    <html>
    <head>
        <title>Neoantigen Report - ${meta.id}</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 40px;
                background-color: #f5f5f5;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background: white;
                padding: 30px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            h1 {{
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
            }}
            h2 {{
                color: #34495e;
                margin-top: 30px;
            }}
            .header-info {{
                background: #ecf0f1;
                padding: 15px;
                border-radius: 5px;
                margin-bottom: 20px;
            }}
            .summary-box {{
                display: inline-block;
                background: #3498db;
                color: white;
                padding: 15px 25px;
                border-radius: 5px;
                margin: 10px 10px 10px 0;
                text-align: center;
            }}
            .summary-box.class-i {{
                background: #2ecc71;
            }}
            .summary-box.class-ii {{
                background: #9b59b6;
            }}
            .summary-number {{
                font-size: 32px;
                font-weight: bold;
            }}
            .summary-label {{
                font-size: 12px;
                opacity: 0.9;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin-top: 15px;
            }}
            th, td {{
                padding: 12px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }}
            th {{
                background-color: #3498db;
                color: white;
            }}
            tr:hover {{
                background-color: #f5f5f5;
            }}
            .hla-box {{
                background: #e8f4f8;
                padding: 15px;
                border-radius: 5px;
                font-family: monospace;
            }}
            .priority-high {{
                color: #27ae60;
                font-weight: bold;
            }}
            .priority-medium {{
                color: #f39c12;
            }}
            .priority-low {{
                color: #95a5a6;
            }}
            .footer {{
                margin-top: 30px;
                padding-top: 20px;
                border-top: 1px solid #ddd;
                color: #7f8c8d;
                font-size: 12px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>🧬 Neoantigen Prediction Report</h1>
            
            <div class="header-info">
                <strong>Sample ID:</strong> ${meta.id}<br>
                <strong>Report Date:</strong> {datetime.now().strftime("%Y-%m-%d %H:%M")}<br>
                <strong>Pipeline:</strong> Neoantigen Pipeline v1.0.0
            </div>
            
            <h2>Summary</h2>
            <div class="summary-box">
                <div class="summary-number">{summary['total_predictions']}</div>
                <div class="summary-label">Total Predictions</div>
            </div>
            <div class="summary-box">
                <div class="summary-number">{summary['filtered_count']}</div>
                <div class="summary-label">Filtered Candidates</div>
            </div>
            <div class="summary-box class-i">
                <div class="summary-number">{summary['class_i_count']}</div>
                <div class="summary-label">Class I</div>
            </div>
            <div class="summary-box class-ii">
                <div class="summary-number">{summary['class_ii_count']}</div>
                <div class="summary-label">Class II</div>
            </div>
            
            <h2>HLA Typing Results</h2>
            <div class="hla-box">
                {hla_alleles.replace(',', ', ')}
            </div>
            
            <h2>Top Neoantigen Candidates</h2>
    '''
    
    if summary['top_candidates']:
        html_content += '''
            <table>
                <tr>
                    <th>Rank</th>
                    <th>Gene</th>
                    <th>Peptide</th>
                    <th>HLA Allele</th>
                    <th>IC50 (nM)</th>
                    <th>%ile Rank</th>
                    <th>Class</th>
                    <th>Score</th>
                </tr>
        '''
        
        for i, candidate in enumerate(summary['top_candidates'], 1):
            score = candidate.get('priority_score', 0)
            if score >= 70:
                score_class = 'priority-high'
            elif score >= 40:
                score_class = 'priority-medium'
            else:
                score_class = 'priority-low'
            
            ic50 = candidate.get('ic50')
            ic50_str = f"{ic50:.1f}" if ic50 else "N/A"
            
            pct = candidate.get('percentile')
            pct_str = f"{pct:.2f}" if pct else "N/A"
            
            html_content += f'''
                <tr>
                    <td>{i}</td>
                    <td><strong>{candidate.get('gene', 'Unknown')}</strong></td>
                    <td style="font-family: monospace;">{candidate.get('peptide', 'Unknown')}</td>
                    <td>{candidate.get('hla_allele', 'Unknown')}</td>
                    <td>{ic50_str}</td>
                    <td>{pct_str}</td>
                    <td>{candidate.get('mhc_class', 'Unknown')}</td>
                    <td class="{score_class}">{score}</td>
                </tr>
            '''
        
        html_content += '</table>'
    else:
        html_content += '<p>No neoantigen candidates passed filtering criteria.</p>'
    
    html_content += '''
            <h2>Methods</h2>
            <p>
                HLA Class I typing was performed using <strong>OptiType</strong> (v1.3.5).
                HLA Class II typing was performed using <strong>HLA-HD</strong> (v1.7.0).
                Neoantigen prediction was performed using <strong>pVACtools</strong> with 
                <strong>NetMHCpan</strong> (Class I) and <strong>NetMHCIIpan</strong> (Class II).
            </p>
            <p>
                Filtering criteria: IC50 ≤ 500 nM, percentile rank ≤ 2%, VAF ≥ 5%.
                Priority scoring considers binding affinity, differential agretopicity, and clonality.
            </p>
            
            <div class="footer">
                Generated by Neoantigen Pipeline | Bioinformatics Core Facility
            </div>
        </div>
    </body>
    </html>
    '''
    
    # Write HTML report
    with open('${meta.id}_neoantigen_report.html', 'w') as f:
        f.write(html_content)
    
    # Log
    with open('${meta.id}_report.log', 'w') as f:
        f.write(f"Report generated for ${meta.id}\\n")
        f.write(f"Top candidates: {len(summary['top_candidates'])}\\n")
    
    print(f"Report generated: ${meta.id}_neoantigen_report.html")
    """
}
