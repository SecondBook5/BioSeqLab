# dna_analyzer/faster_parser.py
# Module for reading and parsing FASTA files

# Import necessary modules from Biopython
from Bio import SeqIO  # Module to help with sequence input/output

class FastaParser:
    def __init__(self, filename):
        """
        Initialize the FastaParser with the input filename.

        Args:
            filename (str): The path to the multi-FASTA file.
        """
        self.filename = filename  # Attribute initialization in an instance variable 'filename' passed to our constructor
        self.records = []  # Initialize an empty list to store sequence records

    def read_sequences(self):
        """
        Read and parse a multi-FASTA file (either in FASTA format or in a FASTQ file) these sequences are stored in records.

        Raises:
            FileNotFoundError: If the file is not found.
            ValueError: If the file has an unsupported extension (anything other than .fasta, .fa, or .txt)
            Exception: For other file reading errors.

        Returns:
            None - since the parsed files are stored in object records, instead of returning the records we have a separate method to access parsed records through get_records method
        """
        self.records = []  # Clear any existing records from our new list

        try:
            if not self.filename.endswith(('.fasta', '.fa', '.txt')):  # Determine that the input file is specifically either a .fasta, .fa, or .txt file
                raise ValueError("Unsupported file extension. Only '.fasta', '.fa', and '.txt' files are supported")  # Raise a value error
            with open(self.filename, 'r') as file:  # Open file and read its contents
                for record in SeqIO.parse(file, "fasta"):  # Iterate through file and parse sequences using Biopython
                    self.records.append({"header": record.id, "sequence": str(record.seq)})  # Append parsed record to list of sequences
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.filename}")  # Handle file not found error
        except Exception as e:
            raise Exception(f"Unexpected Error: {e}")  # Handle other file reading errors

    def get_records(self):
        """
        Get the parsed records.

        Returns:
            list: A list of dictionaries, each containing 'header' and 'sequence'.
        """
        return self.records  # Return the list of parsed records
