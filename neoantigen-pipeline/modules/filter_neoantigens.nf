/*
========================================================================================
    FILTER NEOANTIGENS
========================================================================================
    Apply clinical filtering criteria and rank neoantigen candidates
*/

process FILTER_NEOANTIGENS {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/${meta.id}/neoantigens", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(pvacseq_dir)

    output:
    tuple val(meta), path("${meta.id}_neoantigens_filtered.tsv"),   emit: filtered
    tuple val(meta), path("${meta.id}_neoantigens_ranked.tsv"),     emit: ranked
    tuple val(meta), path("${meta.id}_neoantigen_summary.json"),    emit: summary
    path "*.log",                                                    emit: logs

    script:
    def binding_thresh = params.binding_threshold
    def percentile_thresh = params.percentile_threshold
    def min_vaf = params.min_vaf
    """
    #!/usr/bin/env python3
    
    import os
    import json
    import pandas as pd
    from pathlib import Path
    
    def load_pvacseq_results(pvacseq_dir):
        \"\"\"Load and combine pVACseq results from Class I and II.\"\"\"
        results = []
        
        # Class I results
        class_i_dir = Path(pvacseq_dir) / 'MHC_Class_I'
        if class_i_dir.exists():
            for tsv in class_i_dir.glob('*.all_epitopes.tsv'):
                df = pd.read_csv(tsv, sep='\\t')
                df['MHC_Class'] = 'I'
                results.append(df)
        
        # Class II results
        class_ii_dir = Path(pvacseq_dir) / 'MHC_Class_II'
        if class_ii_dir.exists():
            for tsv in class_ii_dir.glob('*.all_epitopes.tsv'):
                df = pd.read_csv(tsv, sep='\\t')
                df['MHC_Class'] = 'II'
                results.append(df)
        
        if results:
            return pd.concat(results, ignore_index=True)
        return pd.DataFrame()
    
    def calculate_priority_score(row):
        \"\"\"Calculate neoantigen priority score based on multiple factors.\"\"\"
        score = 0
        
        # Binding affinity (lower IC50 = better)
        if 'Median MT IC50 Score' in row and pd.notna(row['Median MT IC50 Score']):
            ic50 = row['Median MT IC50 Score']
            if ic50 < 50:
                score += 40
            elif ic50 < 150:
                score += 30
            elif ic50 < 500:
                score += 20
            elif ic50 < 1000:
                score += 10
        
        # Percentile rank (lower = better)
        if 'Median MT Percentile' in row and pd.notna(row['Median MT Percentile']):
            pct = row['Median MT Percentile']
            if pct < 0.5:
                score += 30
            elif pct < 1:
                score += 20
            elif pct < 2:
                score += 10
        
        # DAI (differential agretopicity index) - MT vs WT binding
        if 'Median WT IC50 Score' in row and 'Median MT IC50 Score' in row:
            if pd.notna(row['Median WT IC50 Score']) and pd.notna(row['Median MT IC50 Score']):
                wt_ic50 = row['Median WT IC50 Score']
                mt_ic50 = row['Median MT IC50 Score']
                if wt_ic50 > 0 and mt_ic50 > 0:
                    dai = wt_ic50 / mt_ic50
                    if dai > 10:
                        score += 20
                    elif dai > 5:
                        score += 15
                    elif dai > 2:
                        score += 10
        
        # VAF (higher = more clonal)
        if 'Tumor DNA VAF' in row and pd.notna(row['Tumor DNA VAF']):
            vaf = row['Tumor DNA VAF']
            if vaf > 0.3:
                score += 10
            elif vaf > 0.1:
                score += 5
        
        return score
    
    # Load results
    df = load_pvacseq_results('${pvacseq_dir}')
    
    if df.empty:
        print("No neoantigen predictions found")
        # Create empty output files
        pd.DataFrame().to_csv('${meta.id}_neoantigens_filtered.tsv', sep='\\t', index=False)
        pd.DataFrame().to_csv('${meta.id}_neoantigens_ranked.tsv', sep='\\t', index=False)
        
        summary = {
            'sample_id': '${meta.id}',
            'total_predictions': 0,
            'filtered_count': 0,
            'class_i_count': 0,
            'class_ii_count': 0,
            'top_candidates': []
        }
        with open('${meta.id}_neoantigen_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        with open('${meta.id}_filter.log', 'w') as f:
            f.write("No neoantigens found for ${meta.id}\\n")
        
        exit(0)
    
    print(f"Loaded {len(df)} total predictions")
    
    # Apply filtering criteria
    filtered = df.copy()
    
    # Binding threshold
    if 'Median MT IC50 Score' in filtered.columns:
        filtered = filtered[filtered['Median MT IC50 Score'] <= ${binding_thresh}]
    
    # Percentile threshold
    if 'Median MT Percentile' in filtered.columns:
        filtered = filtered[filtered['Median MT Percentile'] <= ${percentile_thresh}]
    
    # VAF threshold
    if 'Tumor DNA VAF' in filtered.columns:
        filtered = filtered[filtered['Tumor DNA VAF'] >= ${min_vaf}]
    
    # Remove variants with similar WT binding (not immunogenic)
    if 'Median WT IC50 Score' in filtered.columns and 'Median MT IC50 Score' in filtered.columns:
        # Keep if WT binds poorly OR MT binds much better
        filtered = filtered[
            (filtered['Median WT IC50 Score'] > 500) | 
            (filtered['Median WT IC50 Score'] / filtered['Median MT IC50 Score'] > 1.5)
        ]
    
    print(f"After filtering: {len(filtered)} candidates")
    
    # Calculate priority scores
    filtered['Priority_Score'] = filtered.apply(calculate_priority_score, axis=1)
    
    # Rank by priority score
    ranked = filtered.sort_values('Priority_Score', ascending=False)
    
    # Save filtered results
    filtered.to_csv('${meta.id}_neoantigens_filtered.tsv', sep='\\t', index=False)
    
    # Save ranked results
    ranked.to_csv('${meta.id}_neoantigens_ranked.tsv', sep='\\t', index=False)
    
    # Generate summary
    summary = {
        'sample_id': '${meta.id}',
        'total_predictions': len(df),
        'filtered_count': len(filtered),
        'class_i_count': len(filtered[filtered['MHC_Class'] == 'I']) if 'MHC_Class' in filtered.columns else 0,
        'class_ii_count': len(filtered[filtered['MHC_Class'] == 'II']) if 'MHC_Class' in filtered.columns else 0,
        'filtering_criteria': {
            'binding_threshold': ${binding_thresh},
            'percentile_threshold': ${percentile_thresh},
            'min_vaf': ${min_vaf}
        },
        'top_candidates': []
    }
    
    # Add top 10 candidates to summary
    for _, row in ranked.head(10).iterrows():
        candidate = {
            'gene': row.get('Gene Name', 'Unknown'),
            'mutation': row.get('HGVSc', row.get('Mutation', 'Unknown')),
            'peptide': row.get('MT Epitope Seq', 'Unknown'),
            'hla_allele': row.get('HLA Allele', 'Unknown'),
            'ic50': row.get('Median MT IC50 Score', None),
            'percentile': row.get('Median MT Percentile', None),
            'priority_score': row.get('Priority_Score', 0),
            'mhc_class': row.get('MHC_Class', 'Unknown')
        }
        summary['top_candidates'].append(candidate)
    
    with open('${meta.id}_neoantigen_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Log
    with open('${meta.id}_filter.log', 'w') as f:
        f.write(f"Neoantigen filtering for ${meta.id}\\n")
        f.write(f"Total predictions: {len(df)}\\n")
        f.write(f"After filtering: {len(filtered)}\\n")
        f.write(f"Class I: {summary['class_i_count']}\\n")
        f.write(f"Class II: {summary['class_ii_count']}\\n")
        f.write(f"Top candidate: {summary['top_candidates'][0] if summary['top_candidates'] else 'None'}\\n")
    """
}
