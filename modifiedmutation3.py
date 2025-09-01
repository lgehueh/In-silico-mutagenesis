'''=== enhanced_mutation_generator_deduplicated.py ==='''
import random
from typing import List, Dict, Tuple
import os
import sys

class EnhancedMutationGenerator:
    def __init__(self):
        self.all_mutations = []
        self.codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        self.aa_to_codons = {aa: [] for aa in set(self.codon_table.values())}
        for codon, aa in self.codon_table.items():
            self.aa_to_codons[aa].append(codon)

        self.similar_amino_acids = {
            'E': ['D'], 'D': ['E'], 'K': ['R', 'H'], 'R': ['K', 'H'], 'H': ['K', 'R'],
            'A': ['V', 'G'], 'V': ['A', 'L', 'I'], 'L': ['V', 'I', 'M'],
            'I': ['V', 'L', 'M'], 'M': ['L', 'I', 'V'], 'F': ['Y', 'W'],
            'Y': ['F', 'W'], 'W': ['F', 'Y'], 'S': ['T'], 'T': ['S'],
            'N': ['Q'], 'Q': ['N'], 'C': ['S'], 'G': ['A'], 'P': []
        }
        self.wobble_pairs = {'A': ['T', 'U'], 'T': ['A', 'G'], 'U': ['A', 'G'], 'G': ['C', 'T', 'U'], 'C': ['G']}
        self.wobble_priorities = {2: 1.0, 1: 0.3, 0: 0.1}
        self.near_cognates = self._generate_enhanced_near_cognates()

    def _generate_enhanced_near_cognates(self) -> Dict[str, List[Dict]]:
        near_cognates = {}
        for codon, aa in self.codon_table.items():
            near_cogs = []
            for similar_aa in self.similar_amino_acids.get(aa, []):
                for alt in self.aa_to_codons.get(similar_aa, []):
                    diff = [i for i in range(3) if codon[i] != alt[i]]
                    if len(diff) == 1:
                        i = diff[0]
                        b1, b2 = codon[i], alt[i]
                        wobble = (i == 2) or (b2 in self.wobble_pairs.get(b1, []))
                        if wobble:
                            near_cogs.append({
                                'codon': alt,
                                'amino_acid': similar_aa,
                                'change_position': i,
                                'priority_score': self.wobble_priorities[i] * 1.0,
                                'wobble_mechanism': f'wobble_{i+1}_{b1}->{b2}',
                                'chemical_similarity': 1.0
                            })
            near_cognates[codon] = near_cogs
        return near_cognates

    def get_near_cognate_codons(self, codon: str, include_metadata=False) -> List:
        return self.near_cognates.get(codon, []) if include_metadata else [d['codon'] for d in self.near_cognates.get(codon, [])]

    def read_fasta(self, file: str) -> Dict[str, str]:
        with open(file) as f:
            lines = f.read().splitlines()
        seqs, name, seq = {}, "", ""
        for line in lines:
            if line.startswith('>'):
                if name: seqs[name] = seq
                name, seq = line[1:], ""
            else:
                seq += line.strip().upper()
        if name: seqs[name] = seq
        return seqs

    def parse_sequences(self, sequences: Dict[str, str]) -> Dict[str, Dict[str, str]]:
        return {
            'mrna': {k: v for k, v in sequences.items() if set(v) <= set('ATGC')},
            'protein': {k: v for k, v in sequences.items() if set(v) - set('ATGC')}
        }

    def mrna_pos_to_protein_pos(self, mrna_pos: int) -> int:
        return mrna_pos // 3

    def generate_wobble_mutations(self, mrna_seq, protein_seq, mrna_pos):
        mutations = []
        ref_pos = self.mrna_pos_to_protein_pos(mrna_pos)
        for aa_pos in range(max(0, ref_pos - 21), ref_pos):
            codon = mrna_seq[aa_pos * 3:aa_pos * 3 + 3]
            aa = protein_seq[aa_pos]
            for nc in self.get_near_cognate_codons(codon, include_metadata=True):
                if nc['amino_acid'] != aa and nc['amino_acid'] != '*':
                    new_pep = protein_seq[:aa_pos] + nc['amino_acid'] + protein_seq[aa_pos+1:]
                    mutations.append({
                        'reference_mrna_position': mrna_pos,
                        'reference_protein_position': ref_pos,
                        'mutated_amino_acid_position': aa_pos,
                        'distance_from_reference': ref_pos - aa_pos,
                        'original_codon': codon,
                        'near_cognate_codon': nc['codon'],
                        'original_amino_acid': aa,
                        'near_cognate_amino_acid': nc['amino_acid'],
                        'amino_acid_change': f"{aa}{aa_pos+1}{nc['amino_acid']}",
                        'mutated_peptide': new_pep,
                        'mutation_type': 'enhanced_near_cognate',
                        'change_position': nc['change_position'] + 1,
                        'priority_score': nc['priority_score'],
                        'wobble_mechanism': nc['wobble_mechanism'],
                        'chemical_similarity': nc['chemical_similarity']
                    })
        return mutations

    def generate_near_cognate_mutations(self, mrna_seq, protein_seq, mrna_pos) -> List[Dict]:
        """
        Generate mutations based only on chemical similarity (ignores 3rd position wobble).
        """
        mutations = []
        ref_prot_pos = self.mrna_pos_to_protein_pos(mrna_pos)
        start_pos = max(0, ref_prot_pos - 21)
        end_pos = ref_prot_pos

        for aa_pos in range(start_pos, end_pos):
            codon_start = aa_pos * 3
            if codon_start + 2 >= len(mrna_seq):
                continue

            codon = mrna_seq[codon_start:codon_start+3]
            aa = protein_seq[aa_pos]

            for nc in self.get_near_cognate_codons(codon, include_metadata=True):
                nc_codon = nc['codon']
                nc_aa = nc['amino_acid']

                if nc_aa == aa or nc_aa == '*':
                    continue

                mutated_peptide = protein_seq[:aa_pos] + nc_aa + protein_seq[aa_pos+1:]

                # Filter: only allow chemical similarity-based changes
                if nc['chemical_similarity'] < 1.0:
                    continue

                mutations.append({
                    'reference_mrna_position': mrna_pos,
                    'reference_protein_position': ref_prot_pos,
                    'mutated_amino_acid_position': aa_pos,
                    'distance_from_reference': ref_prot_pos - aa_pos,
                    'original_codon': codon,
                    'near_cognate_codon': nc_codon,
                    'original_amino_acid': aa,
                    'near_cognate_amino_acid': nc_aa,
                    'amino_acid_change': f"{aa}{aa_pos+1}{nc_aa}",
                    'mutated_peptide': mutated_peptide,
                    'mutation_type': 'chemical_similarity_only',
                    'change_position': nc['change_position'] + 1,
                    'priority_score': nc['priority_score'],
                    'wobble_mechanism': nc['wobble_mechanism'],
                    'chemical_similarity': nc['chemical_similarity']
                })
        return mutations

    
    def save_mutations_to_file(self, mutations, path, format_type='fasta', sequence_name="protein"):
        with open(path, 'w') as f:
            if format_type.lower() == 'fasta':
                for i, m in enumerate(mutations):
                    f.write(f">{sequence_name}_mutant_{i+1}_{m['amino_acid_change']}_priority_{m['priority_score']:.2f}_{m['wobble_mechanism']}\n{m['mutated_peptide']}\n")
            else:
                f.write(f"Mutation Results for {sequence_name}\n")
                for i, m in enumerate(mutations):
                    f.write(f"Mutant {i+1}: {m['amino_acid_change']}\n")
                    for key in ['original_codon', 'near_cognate_codon', 'mutation_type', 'wobble_mechanism','original_amino_acid',
                    'near_cognate_amino_acid', 'priority_score', 'chemical_similarity', 'mutated_peptide']:
                        f.write(f"{key}: {m[key]}\n")
                    f.write("-"*40 + "\n")

    def process_sequences(self, fasta_file: str):
        all_sequences = self.read_fasta(fasta_file)
        if not all_sequences:
            print("Error: Could not read input file.")
            return

        parsed = self.parse_sequences(all_sequences)
        if not parsed['mrna'] or not parsed['protein']:
            print("Error: Need both mRNA and protein sequences.")
            return

        mrna_names = list(parsed['mrna'].keys())
        protein_names = list(parsed['protein'].keys())

        print(f"\nAvailable mRNA sequences:")
        for i, name in enumerate(mrna_names, 1):
            print(f"{i}. {name}")

        print(f"\nAvailable protein sequences:")
        for i, name in enumerate(protein_names, 1):
            print(f"{i}. {name}")

        try:
            mrna_input = input("\nSelect mRNA sequence number(s): ")
            protein_input = input("Select protein sequence number(s): ")

            mrna_indices = [int(i.strip()) - 1 for i in mrna_input.split(',')]
            protein_indices = [int(i.strip()) - 1 for i in protein_input.split(',')]

            mutation_mode = input("\nChoose mutation type:\n1. Wobble\n2. Near-cognate\n3. Both\nSelect (1/2/3): ").strip()

            for mrna_idx, protein_idx in zip(mrna_indices, protein_indices):
                mrna_name = mrna_names[mrna_idx]
                protein_name = protein_names[protein_idx]

                mrna_seq = parsed['mrna'][mrna_name]
                protein_seq = parsed['protein'][protein_name]

                print(f"\nSelected mRNA: {mrna_name} (length: {len(mrna_seq)})")
                print(f"Selected Protein: {protein_name} (length: {len(protein_seq)})")

                positions_input = input("\nEnter mRNA position(s) (comma-separated): ")
                mrna_positions = [int(pos.strip()) for pos in positions_input.split(',')]

                os.makedirs("muta-tion", exist_ok=True)

                for mrna_pos in mrna_positions:
                    print(f"\nProcessing position {mrna_pos}...")

                    protein_pos = self.mrna_pos_to_protein_pos(mrna_pos)
                    if protein_pos >= len(protein_seq):
                        print(f"Warning: Protein position {protein_pos} exceeds length.")
                        continue

                    if protein_pos < 21:
                        print(f"Warning: Not enough upstream amino acids before position {protein_pos}")
                        continue

                    self.all_mutations = []

                    if mutation_mode in ['1', '3']:
                        wobble_muts = self.generate_wobble_mutations(mrna_seq, protein_seq, mrna_pos)
                        self.all_mutations.extend(wobble_muts)
                        print(f"Generated {len(wobble_muts)} wobble mutations.")

                    if mutation_mode in ['2', '3']:
                        nc_muts = self.generate_near_cognate_mutations(mrna_seq, protein_seq, mrna_pos)
                        self.all_mutations.extend(nc_muts)
                        print(f"Generated {len(nc_muts)} near-cognate mutations.")

                    seen = set()
                    unique = []
                    for m in self.all_mutations:
                        key = (m['amino_acid_change'], m['mutated_peptide'])
                        if key not in seen:
                            seen.add(key)
                            unique.append(m)

                    if self.all_mutations:
                        base = f"{protein_name}_mRNA{mrna_pos}_protein{protein_pos}"
                        fasta_out = os.path.join("muta-tion", f"{base}.fasta")
                        txt_out = os.path.join("muta-tion", f"{base}.txt")

                        self.save_mutations_to_file(self.all_mutations, fasta_out, format_type='fasta', sequence_name=protein_name)
                        self.save_mutations_to_file(self.all_mutations, txt_out, format_type='txt', sequence_name=protein_name)
                    else:
                        print(f"No mutations generated for position {mrna_pos}.")

        except Exception as e:
            print(f"Error: {e}")

        data = self.read_fasta(fasta_file)
        parsed = self.parse_sequences(data)
        if not parsed['mrna'] or not parsed['protein']:
            print("Error: Missing mRNA or protein sequence.")
            return

        mrna_name = list(parsed['mrna'].keys())[0]
        prot_name = list(parsed['protein'].keys())[0]
        mrna_seq = parsed['mrna'][mrna_name]
        prot_seq = parsed['protein'][prot_name]

        mrna_pos = int(input("Enter mRNA position: "))
        mode = input("Mutation mode (1: Wobble, 2: Near-cognate, 3: Both): ").strip()

        self.all_mutations = []
        if mode in ['1', '3']:
            self.all_mutations.extend(self.generate_wobble_mutations(mrna_seq, prot_seq, mrna_pos))
        if mode in ['2', '3']:
            self.all_mutations.extend(self.generate_wobble_mutations(mrna_seq, prot_seq, mrna_pos))  # Same logic

        # Deduplicate
        seen = set()
        unique = []
        for m in self.all_mutations:
            key = (m['amino_acid_change'], m['mutated_peptide'])
            if key not in seen:
                seen.add(key)
                unique.append(m)

        os.makedirs("deduplicated_output", exist_ok=True)
        base = f"{prot_name}_mRNA{mrna_pos}_protein{self.mrna_pos_to_protein_pos(mrna_pos)}"
        self.save_mutations_to_file(unique, f"deduplicated_output/{base}.fasta", 'fasta', prot_name)
        self.save_mutations_to_file(unique, f"deduplicated_output/{base}.txt", 'txt', prot_name)
        print(f"Saved {len(unique)} unique mutations.")

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 enhanced_mutation_generator_deduplicated.py <input_fasta>")
        return
    EnhancedMutationGenerator().process_sequences(sys.argv[1])

if __name__ == "__main__":
    main()
