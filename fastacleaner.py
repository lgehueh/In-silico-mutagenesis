#!/usr/bin/env python3
"""
FASTA File MaxQuant Cleaner with Systematic Numbering
Processes FASTA files to create systematic numbering for protein variants:
- sp|P24394_1|IL4R_HUMAN_1 ... SV=1_original
- sp|P24394_2|IL4R_HUMAN_2 ... SV=2_P_598_S_1st  
- sp|P24394_3|IL4R_HUMAN_3 ... SV=3_M_120_I_wobble2
- Then moves to next protein: sp|KMT2B_1|KMT2B_HUMAN_1 ... SV=1_original
"""

import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict

def extract_mutation_info(header):
    """
    Extract mutation and wobble information from header
    """
    mutation_info = ""
    wobble_info = ""
    
    # Look for mutation patterns like _P_598_S_1st, _M_120_I, etc.
    mutation_match = re.search(r'(_[A-Z]_\d+_[A-Z](?:_\d+[a-z]*)?)', header)
    if mutation_match:
        mutation_info = mutation_match.group(1).lstrip('_')
    
    # Look for _mut patterns
    mut_match = re.search(r'(_mut\d*)', header)
    if mut_match:
        if mutation_info:
            mutation_info += mut_match.group(1)
        else:
            mutation_info = mut_match.group(1).lstrip('_')
    
    # Look for wobble info (1, 2, or 3 - usually at the end)
    wobble_match = re.search(r'_([123])(?:\s|$)', header)
    if wobble_match:
        wobble_info = f"_wobble{wobble_match.group(1)}"
    
    return mutation_info, wobble_info

def extract_protein_info(header):
    """
    Extract protein information from header
    """
    # Remove the leading '>'
    header = header.lstrip('>')
    
    # Extract UniProt ID 
    uniprot_match = re.match(r'sp\|([^|]+)\|', header)
    if uniprot_match:
        uniprot_id = uniprot_match.group(1)
    else:
        # Fallback for other formats
        parts = header.split()
        if parts:
            first_part = parts[0]
            if '|' in first_part:
                uniprot_id = first_part.split('|')[1] if len(first_part.split('|')) > 1 else first_part
            else:
                uniprot_id = first_part
        else:
            uniprot_id = "UNKNOWN"
    
    # Extract description and metadata
    parts = header.split()
    description_parts = []
    metadata_parts = []
    
    # Look for where metadata starts (OS=, OX=, GN=, PE=, SV=)
    metadata_started = False
    
    for part in parts[1:]:
        if re.match(r'^(OS|OX|GN|PE|SV)=', part):
            metadata_started = True
        
        if metadata_started:
            metadata_parts.append(part)
        else:
            description_parts.append(part)
    
    # Clean description - keep only the first meaningful part
    clean_description = []
    if description_parts:
        for i, part in enumerate(description_parts):
            clean_description.append(part)
            # Stop at reasonable description length or second protein name indicators
            if i > 4 or (i > 1 and (part.isupper() or 
                         any(indicator in part.lower() for indicator in ['variant', 'isoform', 'splice']))):
                break
    
    description = ' '.join(clean_description) if clean_description else "Protein"
    
    # Parse metadata
    metadata = {}
    for part in metadata_parts:
        if '=' in part:
            key, value = part.split('=', 1)
            metadata[key] = value
    
    # Use GN= field as the gene name (this should be consistent)
    gene_name = metadata.get('GN', 'UNKNOWN')
    
    return uniprot_id, gene_name, description, metadata

def process_fasta_file(input_file, output_file=None):
    """
    Process a FASTA file with systematic numbering
    """
    if output_file is None:
        output_file = input_file.with_suffix('.maxquant.fasta')
    
    # Store sequences grouped by protein
    protein_sequences = defaultdict(list)
    
    # First pass: collect all sequences and group by protein
    with open(input_file, 'r') as infile:
        current_header = None
        current_sequence = []
        
        for line in infile:
            line = line.strip()
            
            if line.startswith('>'):
                # Process previous sequence if exists
                if current_header is not None:
                    uniprot_id, gene_name, description, metadata = extract_protein_info(current_header)
                    mutation_info, wobble_info = extract_mutation_info(current_header)
                    
                    # Use gene name as the primary grouping key
                    protein_key = gene_name
                    
                    protein_sequences[protein_key].append({
                        'original_header': current_header,
                        'uniprot_id': uniprot_id,
                        'gene_name': gene_name,
                        'description': description,
                        'metadata': metadata,
                        'mutation_info': mutation_info,
                        'wobble_info': wobble_info,
                        'sequence': current_sequence
                    })
                
                # Start new sequence
                current_header = line
                current_sequence = []
                
            else:
                # Add to current sequence
                current_sequence.append(line)
        
        # Process the last sequence
        if current_header is not None:
            uniprot_id, gene_name, description, metadata = extract_protein_info(current_header)
            mutation_info, wobble_info = extract_mutation_info(current_header)
            
            protein_key = gene_name
            
            protein_sequences[protein_key].append({
                'original_header': current_header,
                'uniprot_id': uniprot_id,
                'gene_name': gene_name,
                'description': description,
                'metadata': metadata,
                'mutation_info': mutation_info,
                'wobble_info': wobble_info,
                'sequence': current_sequence
            })
    
    # Second pass: write with systematic numbering
    with open(output_file, 'w') as outfile:
        for protein_name in sorted(protein_sequences.keys()):
            sequences = protein_sequences[protein_name]
            
            # Sort sequences: original first, then mutations, then wobbles
            def sort_key(seq):
                mutation = seq['mutation_info']
                wobble = seq['wobble_info']
                
                if not mutation and not wobble:
                    return (0, '')  # Original first
                elif mutation and not wobble:
                    return (1, mutation)  # Mutations second
                else:
                    return (2, mutation + wobble)  # Wobbles last
            
            sequences.sort(key=sort_key)
            
            # Write sequences with systematic numbering
            for i, seq in enumerate(sequences, 1):
                # Create systematic header
                new_uniprot_id = f"{seq['uniprot_id']}_{i}"
                new_gene_name = f"{seq['gene_name']}_{i}"
                
                # Create SV field with mutation info
                if seq['mutation_info'] or seq['wobble_info']:
                    sv_value = f"{i}_{seq['mutation_info']}{seq['wobble_info']}"
                else:
                    sv_value = f"{i}_original"
                
                # Build metadata
                metadata = seq['metadata'].copy()
                metadata['SV'] = sv_value
                
                # Construct new header
                header_parts = [f"sp|{new_uniprot_id}|{new_gene_name}", seq['description']]
                
                # Add metadata in standard order
                for key in ['OS', 'OX', 'GN', 'PE', 'SV']:
                    if key in metadata:
                        header_parts.append(f"{key}={metadata[key]}")
                
                new_header = ' '.join(header_parts)
                
                # Write sequence
                outfile.write(f">{new_header}\n")
                for seq_line in seq['sequence']:
                    outfile.write(f"{seq_line}\n")
    
    print(f"Processed {input_file}")
    print(f"Output written to {output_file}")
    print(f"Processed {len(protein_sequences)} different proteins")
    total_sequences = sum(len(sequences) for sequences in protein_sequences.values())
    print(f"Total sequences: {total_sequences}")

def main():
    parser = argparse.ArgumentParser(description='Clean FASTA files with systematic numbering for MaxQuant')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output file (default: input_file.maxquant.fasta)')
    
    args = parser.parse_args()
    
    input_path = Path(args.input_file)
    output_path = Path(args.output) if args.output else None
    
    if not input_path.exists():
        print(f"Error: Input file {input_path} does not exist")
        sys.exit(1)
    
    try:
        process_fasta_file(input_path, output_path)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
