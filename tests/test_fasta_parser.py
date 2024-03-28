#test_fasta_parser.py
#contains test for the FASTA file parser
import unittest
from dna_analyzer.fasta_parser import FastaParser 

class TestFastaParser(unittest.TestCase):
    def setUp(self):
        # Define a common file path for testing
        self.fasta_file = 'path/to/test.fasta'
        self.fastq_file = 'path/to/test.fastq'

    def test_read_sequences_fasta(self):
        parser = FastaParser(self.fasta_file)
        parser.read_sequences()
        records = parser.get_records()
        self.assertTrue(records)  # Check if records are not empty

    def test_read_sequences_fastq(self):
        parser = FastaParser(self.fastq_file)
        parser.read_sequences()
        records = parser.get_records()
        qualities = parser.get_qualities()
        self.assertTrue(records)  # Check if records are not empty
        self.assertTrue(qualities)  # Check if qualities are not empty for FASTQ

    def test_calculate_gc_content(self):
        parser = FastaParser(self.fasta_file)
        parser.read_sequences()
        parser.calculate_gc_content()
        records = parser.get_records()
        for record in records:
            self.assertIn('gc_content', record)  # Check if gc_content is calculated for each record

    def test_extract_subsequence(self):
        parser = FastaParser(self.fasta_file)
        parser.read_sequences()
        parser.extract_subsequence(1, 5)
        records = parser.get_records()
        for record in records:
            self.assertIn('subsequence', record)  # Check if subsequence is extracted for each record

    def test_search_sequences(self):
        parser = FastaParser(self.fasta_file)
        parser.read_sequences()
        keyword = 'ATG'
        found_sequences = parser.search_sequences(keyword)
        self.assertTrue(found_sequences)  # Check if sequences containing the keyword are found

if __name__ == '__main__':
    unittest.main()
