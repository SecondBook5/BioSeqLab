#dna_analyzer/faster_parser.py
#Module for reading and parsing FASTA files

#import necessary modules from Biopython 
from Bio import SeqIO #module to help with sequence input/output
from Bio.Seq import Seq #module helpful for working with biological sequences 
from Bio.SeqRecord import SeqRecord #module useful for creating a container for sequences with metadata


class FastaParser:
    def __init__(self, filename): 
        """
        Initialize the FastaParser with the input filename.

        Args:
            filename (str): The path to the multi-FASTA file.
        """
        self.filename = filename #attribute initialization in an instance variable 'filename' passed to our constructor
        self.records = []  #Initialize an empty list to store sequence records               
        self.qualities = [] #Initialize an empty list to stor qualities if the file is of FASTQ format
    
    def read_sequences(self):
        """
        Read and parse a multi-FASTA file (either in FASTA format or in a FASTQ file) these sequences are stored in records.

        Raises:
            FileNotFoundError: If the file is not found.
            ValueError: If the file has an unsupported extension (anything other than .fasta or .fa)
            Exception: For other file reading errors.

        Returns:
            None- since the parsed files are stored in object records, instead of returning the records we have a separate method to access parsed records through get_records method
        """ 
    
        self.records = []   #Clear any existing records from our new list
        self.qualities = [] #Clear any existing qualities from our new list
        
        # Determine the file extension if its a FASTA file 
        try:
            if not self.filename.endswith(('.fasta', '.fa')): #determine that the input file is specifically either a .fasta or a .fa file
                raise ValueError("Unsupported file extension. Only '.fasta' and '.fa' files are supported") #else raise a value error
            with open(self.filename,'r') as file: #open FASTA file and read its contents
                for record in SeqIO.parse(file,"fasta"): #loop to iterate through FASTA file and parse through the sequences within the FASTA file using Biopython
                    self.record_count.append({"header": record.id, "sequence": str(record.seq)}) #append each parsed record to the list of sequences created in record_count
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.filename}") #handles a file not found error 
        except Exception as e:
            raise Exception(f"Unexpected Error:{e}") #handle any other file error exceptions
    def get_records(self):
        """
        Get the parsed records.

        Returns:
            list: A list of dictionaries, each containing 'header' and 'sequence'.
        """
        return self.record_count #returns the list of parsed records
    






