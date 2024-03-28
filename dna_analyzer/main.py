#main.py
#Entry point code for the program
import argparse
from Bio import SeqIO #Biopython for FASTA parsing
from dna_analyzer.fasta_parser import FastaParser

def main():
    #create an ArgumentParser instance that defines a program description
    parser = argparse.ArgumentParser(description = "DNA Sequence Analysis Tool")

    #Define the parameters of the command-line arguments
    parser.add_argument("input_file", type=str, help="Path to the input multi-FASTA file") #file input command line
    
    
    #Parse the command-line arguments
    args = parser.parse_args()

    #Read and parse the input multi-FASTA file by creating an instance of FastaParser
    fasta_parser = FastaParser() #create an object named fasta_parser from the FastaParser class. This object will parse your multi-FASTA files
    fasta_parser.read_fasta(args.input_file) #uses the read_fasta function from fasta_parser to read the input file entered in path these are then stored in fasta_parser object
    records = fasta_parser.get_records() #this will then return the parsed records 



'''

from orf_finder import ORFFinder
from repeat_finder import RepeatFinder


    parser.add_argument("--frame", type = int, help = "Reading frame for ORF analysis accepts frames for (1, 2, or 3)") #frame command line
    parser.add_argument("--repeat_length", type = int, help = "Length of repeats for analysis") #repeat length command line
    parser.add_argument("--sequence_identifier", type = str, help = "Identifier of a specific sequence") #specific sequence in command line

    

   

    #Create instances of ORFFinder and RepeatFinder: create 2 objects orf_finder and repeat_finder from their respective classes
    #These will then pass through those objects the records list as an argument thus reading the parsed DNA sequences
    orf_finder = ORFFinder(records) 
    repeat_finder = RepeatFinder(records)

    if args.sequence_identifier:
        longest_orf_length_sequence = orf_finder.find_longest_orf_in_sequence(args.sequence_identifier)
        print(f"")
'''