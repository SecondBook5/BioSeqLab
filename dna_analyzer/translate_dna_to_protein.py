from Bio.Seq import Seq

# Define the DNA sequence
dna_sequence = "TACCGTTAACCT"

# Create a Seq object
dna_seq = Seq(dna_sequence)

# Translate the DNA sequence into amino acids (ORF 1)
orf1_protein_seq = dna_seq.translate()

# Translate the reverse complement of the DNA sequence into amino acids (ORF 2)
orf2_protein_seq = dna_seq.reverse_complement().translate()

# Translate the second reading frame into amino acids (ORF 3)
orf3_protein_seq = dna_seq[1:].translate()

# Translate the reverse complement of the second reading frame into amino acids (ORF 4)
orf4_protein_seq = dna_seq.reverse_complement()[1:].translate()

# Translate the third reading frame into amino acids (ORF 5)
orf5_protein_seq = dna_seq[2:].translate()

# Translate the reverse complement of the third reading frame into amino acids (ORF 6)
orf6_protein_seq = dna_seq.reverse_complement()[2:].translate()

# Print the results
print("ORF 1 Protein Sequence:", orf1_protein_seq)
print("ORF 2 Protein Sequence:", orf2_protein_seq)
print("ORF 3 Protein Sequence:", orf3_protein_seq)
print("ORF 4 Protein Sequence:", orf4_protein_seq)
print("ORF 5 Protein Sequence:", orf5_protein_seq)
print("ORF 6 Protein Sequence:", orf6_protein_seq)
