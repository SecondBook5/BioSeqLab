#orf_finder.py
#module for finding Open Reading Frames (ORFs) containing functions for finding the longest ORFs (Open Reading Frames)
"""
https://en.wikipedia.org/wiki/Open_reading_frame 
Definition: an ORF is a sequence that has a length divisible by three and begins with a translation start codon (ATG) and ends at a stop codon. Diff
Further Reading:  https://www-sciencedirect-com.ezproxy.yu.edu/science/article/pii/S0168952517302299


"""
# orf_finder.py
# Module for finding Open Reading Frames (ORFs)
from orfipy import findorfs
from Bio.Seq import Seq, SeqRecord, SeqIO  # Import Seq from Biopython for working with sequences

class ORFFinder:
    def __init(self):
        self.records = []  # Initialize an empty list to store parsed records.

    def set_records(self, records):
        """
        Set the parsed records for ORF analysis.

        Args:
            records (list): A list of dictionaries, each containing 'header' and 'sequence'.
        """
        self.records = records  # Store the parsed records for ORF analysis.

    def find_longest_orf(self, reading_frame):
        """
        Find the length of the longest ORF in the specified reading frame.

        Args:
            reading_frame (int): The reading frame (1, 2, or 3).

        Returns:
            tuple: A tuple containing the length, identifier, and starting position of the longest ORF.
        """
        #Error handling to ensure the reading frame is either 1, 2, 3 or -1, -2, -3 else it returns a value error 
        if reading_frame not in (-3, -2, -1, 1, 2, 3): 
            raise ValueError("Invalid reading frame. The  forward reading frame must be either 1, 2, or 3 and the reverse reading frame must be either -1, -2, -3") 
        
        longest_length = 0  # Initialize a variable to track the length of the longest ORF.
        longest_identifier = None  # Initialize a variable to store the identifier of the longest ORF.
        longest_start = None  # Initialize a variable to store the starting position of the longest ORF.
    
    
    #If reading_frame < 0 to check if reading_frame is negative. If it's negative, it calculates the reverse complement of the sequence and then looks for ORFs. 
    #If reading_frame is positive, it directly looks for ORFs in the original sequence.
        
        for record in self.records:  # Loop through the parsed records.
            sequence = Seq(record['sequence'])  # Convert the DNA sequence to a Seq object.
            #Determine if it is the forward reading frame or the reverse reading frame
            if reading_frame < 0: 
                sequence = sequence.reverse_complement() #if it is the reverse reading frame use reverse complement to determine the reverse ORF
            for orf in self._find_orfs(sequence, abs(reading_frame)):  #begin a loops that iterates over a collection of ORFs we created in the later method _find_orfs. It iterates over each orf in this collection
                if len(orf) > longest_length:
                    longest_length = len(orf)
                    longest_identifier = record['header']
                    longest_start = orf.location.start + 1
        return longest_length, longest_identifier, longest_start
                
    def find_longest_orf_in_file(self):
        """
        Find the length of the longest ORF in any sequence and any forward reading frame.

        Returns:
            tuple: A tuple containing the length, identifier, and starting position of the longest ORF.
        """
        longest_length = 0  # Initialize a variable to track the length of the longest ORF.
        longest_identifier = None  # Initialize a variable to store the identifier of the longest ORF.
        longest_start = None  # Initialize a variable to store the starting position of the longest ORF.
        longest_frame = None  # Initialize a variable to store the reading frame of the longest ORF.

        for record in self.records:  # Loop through the parsed records.
            for reading_frame in range(1, 4):  # Iterate through reading frames (1, 2, 3).
                orf_length, _, _ = self.find_longest_orf(reading_frame)  # Find the longest ORF in the current reading frame.
                if orf_length > longest_length:  # Compare with the longest found so far.
                    longest_length = orf_length  # Update the length of the longest ORF.
                    longest_identifier = record['header']  # Update the identifier of the longest ORF.
                    longest_start = None  # Update this with the actual start position.
                    longest_frame = reading_frame  # Update the reading frame of the longest ORF.

        return longest_length, longest_identifier, longest_start, longest_frame  # Return information about the longest ORF.

    def find_longest_orf_in_sequence(self, sequence_identifier, reading_frame):
        """
        Find the length of the longest forward ORF in a specific sequence.

        Args:
            sequence_identifier (str): The identifier of the sequence to analyze.

        Returns:
            tuple: A tuple containing the length and starting position of the longest ORF.
        """
        if reading_frame not in [-3, -2, -1, 1, 2, 3]:
            raise ValueError("Invalid reading frame. The  forward reading frame must be either 1, 2, or 3 and the reverse reading frame must be either -1, -2, -3")
        
        longest_length = 0  # Initialize a variable to track the length of the longest ORF.
        longest_start = None  # Initialize a variable to store the starting position of the longest ORF.

        for record in self.records:  # Loop through the parsed records.
            if record['header'] == sequence_identifier:  # Check if the record matches the specified identifier.
                sequence = Seq.Seq(record['sequence'])  # Convert the DNA sequence to a Seq object.
            for orf in self._find_orfs(sequence, abs(reading_frame)):
                if len(orf) > longest_length:
                        longest_length = len(orf)
                        longest_start = orf.location.start + 1

        return longest_length, longest_start
                
    
    def _find_orfs(self, sequence, reading_frame):
        """
        Find and collect valid ORFs (Open Reading Frames) in a given DNA sequence.

        Args:
            sequence (Seq): The DNA sequence to search for ORFs.
            reading_frame (int): The reading frame (1, 2, 3 for forward, or -1, -2, -3 for reverse).

        Returns:
            list: A list of valid ORFs found in the specified reading frame.
        """
        # Create an empty list to store identified ORFs.
        orf_list = []

        # Iterate through the DNA sequence, searching for start codons (ATG).
        for orf in SeqIO.SeqUtils.nt_search(sequence, Seq.Seq("ATG"), start=reading_frame - 1):
            # Translate the found ORF into an amino acid sequence.
            translation = orf.translate()

            # Check if the translated sequence contains a stop codon ('*') at the end.
            if '*' not in translation[:-1]:
                # If there's no stop codon, it's a valid ORF. Append it to the list.
                orf_list.append(orf)

        # Return the list of valid ORFs.
        return orf_list






'''
            start_range = reading_frame -1 #forward reading frame starts in position 0, 1 or 2
                step = 3 #count 3 codons at a time moving forward
            else: 
                start_range = len(sequence) + reading_frame #length of the sequence + whatever orf length you are starting off with
                step = -3

            for i in range(start_range, len(sequence), step):  # Iterate through the sequence based on the reading frame.
                codon = str(sequence[i:i + 3])  # Extract a codon of three nucleotides.
                if codon == 'ATG':  # Check if the codon is a start codon (ATG).
                    orf_length = 0  # Initialize a variable to track the ORF length.
                    for j in range(i, len(sequence), step):  # Continue iterating through the sequence.
                        codon = str(sequence[j:j + 3])  # Extract the next codon.
                        if codon in ('TAA', 'TAG', 'TGA'):  # Check if the codon is a stop codon.
                            break  # Exit the loop if a stop codon is found.
                        orf_length += 3  # Increment the ORF length by 3 for each codon.
                    if orf_length > longest_length:  # Compare the length of the current ORF with the longest found so far.
                        longest_length = orf_length  # Update the length of the longest ORF.
                        longest_identifier = record['header']  # Update the identifier of the longest ORF.
                        longest_start = i + 1  # Update the starting position of the longest ORF.
        return longest_length, longest_identifier, longest_start  # Return information about the longest ORF.
'''
'''
for reading_frame in range(1, 4):  # Iterate through reading frames (1, 2, 3).
                    orf_length, start, _ = self.find_longest_orf(reading_frame)  # Find the longest ORF in the current reading frame.
                    if orf_length > longest_length:  # Compare with the longest found so far.
                        longest_length = orf_length  # Update the length of the longest ORF.
                        longest_start = start  # Update the starting position of the longest ORF.

        return longest_length, longest_start  # Return information about the longest ORF in the specified sequence.
'''