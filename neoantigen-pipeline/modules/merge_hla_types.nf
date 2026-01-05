/*
========================================================================================
    MERGE HLA TYPES
========================================================================================
    Combine HLA Class I (OptiType) and Class II (HLA-HD) results into unified format
*/

process MERGE_HLA_TYPES {
    tag "$meta.id"
    label 'process_low'
    
    publishDir "${params.outdir}/${meta.id}/hla_typing", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(optitype_result), path(hlahd_result)

    output:
    tuple val(meta), path("${meta.id}_hla_alleles.txt"),    emit: merged_hla
    tuple val(meta), path("${meta.id}_hla_summary.json"),   emit: hla_json
    path "*.log",                                            emit: logs

    script:
    """
    #!/usr/bin/env python3
    
    import json
    import re
    import sys
    
    def parse_optitype(filepath):
        \"\"\"Parse OptiType TSV output for Class I alleles.\"\"\"
        alleles = {'A': [], 'B': [], 'C': []}
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('\\t') or line.startswith('0\\t'):
                    parts = line.strip().split('\\t')
                    # OptiType format: idx, A1, A2, B1, B2, C1, C2, Reads, Objective
                    if len(parts) >= 7:
                        alleles['A'] = [parts[1], parts[2]]
                        alleles['B'] = [parts[3], parts[4]]
                        alleles['C'] = [parts[5], parts[6]]
        return alleles
    
    def parse_hlahd(filepath):
        \"\"\"Parse HLA-HD result file for all HLA alleles.\"\"\"
        alleles = {}
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\\t')
                if len(parts) >= 2:
                    gene = parts[0].replace('HLA-', '')
                    # HLA-HD reports two alleles per gene
                    typed_alleles = [a for a in parts[1:] if a and a != '-' and a != 'Not typed']
                    if typed_alleles:
                        alleles[gene] = typed_alleles[:2]  # Max 2 alleles
        return alleles
    
    def format_allele_for_pvacseq(allele):
        \"\"\"Format HLA allele for pVACseq (HLA-A*02:01 format).\"\"\"
        # Remove 'HLA-' prefix if present, then add it back in standard format
        allele = allele.replace('HLA-', '')
        # Ensure proper format: A*02:01:01 -> HLA-A*02:01
        match = re.match(r'([A-Z]+[0-9]?)\\*?(\\d+):(\\d+)', allele)
        if match:
            gene, group, protein = match.groups()
            return f"HLA-{gene}*{group}:{protein}"
        return None
    
    # Parse results
    class_i = parse_optitype('${optitype_result}')
    class_ii = parse_hlahd('${hlahd_result}')
    
    # Use OptiType for Class I (more accurate for A, B, C)
    # Use HLA-HD for Class II (DRB1, DQB1, DPB1, etc.)
    
    all_alleles = {
        'class_i': {
            'A': class_i.get('A', []),
            'B': class_i.get('B', []),
            'C': class_i.get('C', [])
        },
        'class_ii': {
            'DRB1': class_ii.get('DRB1', []),
            'DQA1': class_ii.get('DQA1', []),
            'DQB1': class_ii.get('DQB1', []),
            'DPA1': class_ii.get('DPA1', []),
            'DPB1': class_ii.get('DPB1', [])
        }
    }
    
    # Format alleles for pVACseq
    pvacseq_alleles = []
    
    # Class I alleles
    for gene in ['A', 'B', 'C']:
        for allele in all_alleles['class_i'].get(gene, []):
            formatted = format_allele_for_pvacseq(allele)
            if formatted and formatted not in pvacseq_alleles:
                pvacseq_alleles.append(formatted)
    
    # Class II alleles
    for gene in ['DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']:
        for allele in all_alleles['class_ii'].get(gene, []):
            formatted = format_allele_for_pvacseq(allele)
            if formatted and formatted not in pvacseq_alleles:
                pvacseq_alleles.append(formatted)
    
    # Write alleles file for pVACseq (comma-separated)
    with open('${meta.id}_hla_alleles.txt', 'w') as f:
        f.write(','.join(pvacseq_alleles))
    
    # Write JSON summary
    summary = {
        'sample_id': '${meta.id}',
        'class_i_source': 'OptiType',
        'class_ii_source': 'HLA-HD',
        'alleles': all_alleles,
        'pvacseq_format': pvacseq_alleles
    }
    
    with open('${meta.id}_hla_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    # Log
    with open('${meta.id}_merge_hla.log', 'w') as f:
        f.write(f"Merged HLA types for ${meta.id}\\n")
        f.write(f"Class I alleles: {all_alleles['class_i']}\\n")
        f.write(f"Class II alleles: {all_alleles['class_ii']}\\n")
        f.write(f"pVACseq format: {','.join(pvacseq_alleles)}\\n")
    
    print(f"Found {len(pvacseq_alleles)} HLA alleles for ${meta.id}")
    """
}
