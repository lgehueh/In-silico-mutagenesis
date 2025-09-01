# In-silico-mutagenesis
This workflow generates mutant variants of the candidate canonical protein sequences.

gen_mutation.py – Accepts a FASTA file of protein sequences as input. The user specifies a G4 position, and the script introduces mutations within a 21-amino-acid window upstream of that position. Mutations are generated according to first, second, and third wobble logic, and the resulting mutant sequences are written to a new FASTA file.

fastacleaner.py – Takes FASTA files as input, removes redundant sequences, and organizes the protein entries with unique IDs in a consistent order. This ensures compatibility with MaxQuant, which requires all sequence entries to have non-redundant identifiers.
