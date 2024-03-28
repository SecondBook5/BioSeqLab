#repeat_finder.py
#module includes the function for finding DNA repeats
from Bio.Seq import Seq

class RepeatFinder:
    def __init__(self):
        self.records = []  # Initialize an empty list to store parsed records.

    def set_records(self, records):
        """
        Set the parsed records for repeat analysis.

        Args:
            records (list): A list of dictionaries, each containing 'header' and 'sequence'.
        """
        self.records = records  # Store the parsed records for repeat analysis.

    def find_repeats(self, length):
        """
        Find repeats of a specific length in all sequences.

        Args:
            length (int): The length of repeats to search for.

        Returns:
            list: A list of dictionaries, each containing 'sequence' and 'count'.
        """
        repeats = {}  # Initialize a dictionary to store repeats and their counts.

        for record in self.records:  # Loop through the parsed records.
            sequence = Seq(record['sequence'])  # Convert the DNA sequence to a BioPython Seq object. We use BioPython's `Seq` class to represent DNA sequences, making it easier to work with sequences and substrings.
            sequence_str = str(sequence)  # Convert the Seq object to a string for substring extraction.
            for i in range(len(sequence_str) - length + 1):  # Iterate through the sequence.
                repeat = sequence_str[i:i + length]  # Extract a substring of the specified length.
                if repeat in repeats:
                    repeats[repeat] += 1  # Increment the repeat count.
                else:
                    repeats[repeat] = 1  # Initialize the repeat count.

        repeat_list = [{'sequence': repeat, 'count': count} for repeat, count in repeats.items()]  # Convert repeats to a list of dictionaries.

        return repeat_list  # Return a list of repeats and their counts.

    def find_most_frequent_repeat(self, length):
        """
        Find the most frequently occurring repeat of a specific length in all sequences.

        Args:
            length (int): The length of repeats to search for.

        Returns:
            tuple: A tuple containing the sequence and count of the most frequent repeat.
        """
        repeats = self.find_repeats(length)  # Find all repeats of the specified length.
        most_frequent_repeat = max(repeats, key=lambda x: x['count'])  # Find the most frequent repeat.

        return most_frequent_repeat  # Return the most frequent repeat.

    def find_most_frequent_repeats(self, length, max_count):
        """
        Find the most frequently occurring repeats of a specific length in all sequences.

        Args:
            length (int): The length of repeats to search for.
            max_count (int): The maximum number of most frequent repeats to find.

        Returns:
            list: A list of dictionaries, each containing 'sequence' and 'count' of the most frequent repeats.
        """
        repeats = self.find_repeats(length)  # Find all repeats of the specified length.
        most_frequent_repeats = sorted(repeats, key=lambda x: x['count'], reverse=True)[:max_count]  # Find the most frequent repeats.

        return most_frequent_repeats  # Return a list of the most frequent repeats.
    
