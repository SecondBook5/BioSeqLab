#
from Bio import SeqIO
import logging
from bioseqlab.parsers.sequence_parser import SequenceParser, SequenceParserError


class FastaParser(SequenceParser):
    """
    FastaParser is a class for parsing FASTA sequence files.

    This class provides methods to read sequences from a FASTA file,
    retrieve headers and sequences, plot sequence length distributions,
    and compute sequence statistics.
    """

    @property
    def supported_extensions(self) -> list:
        """
        Get the list of supported file extensions for FASTA files.

        Returns:
            list: A list of supported file extensions, such as ['.fasta', '.fa'].
        """
        # Return the list of supported file extensions for FASTA files
        return ['.fasta', '.fa']

    @staticmethod
    def _validate_nucleotide_sequence(sequence: str) -> bool:
        """
        Validate that the sequence contains only valid nucleotide characters (A, T, C, G, N).

        Args:
            sequence (str): The sequence to validate.

        Returns:
            bool: True if the sequence is valid, False otherwise.
        """
        # Define the set of valid nucleotide characters
        valid_nucleotides = {'A', 'T', 'C', 'G', 'N'}
        # Check if all characters in the sequence are valid nucleotides
        return all(nuc in valid_nucleotides for nuc in sequence)

    def read_sequences(self) -> None:
        """
        Read sequences from the FASTA file.

        Raises:
            SequenceParserError: If an error occurs while reading the sequences.
        """
        # Validate the file type to ensure it has a supported extension
        self._validate_file_type(self.supported_extensions)

        try:
            # Open the FASTA file for reading
            with open(self.filename, 'r') as file:
                self.sequences = []  # Initialize an empty list to store sequences
                # Parse each record in the FASTA file
                for record in SeqIO.parse(file, "fasta"):
                    sequence = str(record.seq)  # Convert the sequence to a string
                    # Validate the nucleotide sequence
                    if not self._validate_nucleotide_sequence(sequence):
                        # Log a warning if the sequence contains invalid characters
                        logging.warning(f"Invalid characters found in sequence with header: {record.id}")
                        continue  # Skip invalid sequences
                    # Add the valid sequence to the sequences list
                    self.sequences.append({"header": record.id, "sequence": sequence})
                    # Logging for unusually short sequences
                    if len(sequence) < 50:
                        logging.warning(f"Sequence {record.id} is unusually short: {len(sequence)} bases")
                    # Logging for unusually long sequences
                    elif len(sequence) > 5000:
                        logging.warning(f"Sequence {record.id} is unusually long: {len(sequence)} bases")

            # Check if no valid sequences were found
            if not self.sequences:
                logging.error(f"No valid sequences found in {self.filename}.")
                raise SequenceParserError(f"No valid sequences found in {self.filename}.")

            # Validate the parsed sequences to ensure they meet expected criteria
            self._validate_sequences()
            # Log the success message with the number of sequences parsed
            logging.info(f"Successfully parsed {len(self.sequences)} sequences from {self.filename}.")
        except Exception as e:
            # Log an error message if there was an issue reading the sequences
            logging.error(f"Error reading sequences from {self.filename}: {e}")
            # Raise a custom SequenceParserError with the original exception
            raise SequenceParserError(f"Error reading sequences from {self.filename}: {e}") from e

    def get_headers(self) -> list:
        """
        Get the headers of the sequences.

        This method returns a list of headers (identifiers) for the sequences in the FASTA file.

        Returns:
            list: A list of sequence headers.
        """
        # Extract and return the list of headers from the parsed sequences
        return [seq["header"] for seq in self.sequences]

    def get_sequence_by_header(self, header: str) -> str:
        """
        Get the sequence corresponding to a given header.

        Args:
            header (str): The header of the sequence.

        Returns:
            str: The sequence corresponding to the header.

        Raises:
            SequenceParserError: If the header is not found.
        """
        # Iterate through the sequences to find the one with the matching header
        for seq in self.sequences:
            if seq["header"] == header:
                # Return the sequence if the header matches
                return seq["sequence"]
        # Log an error if the header was not found
        logging.error(f"Header '{header}' not found in sequences.")
        # Raise a custom SequenceParserError if the header was not found
        raise SequenceParserError(f"Header '{header}' not found in sequences.")

    def filter_sequences_by_length(self, min_length: int = 0, max_length: int = None) -> list:
        """
        Filter sequences by minimum and/or maximum length.

        Args:
            min_length (int): Minimum length of sequences to include.
            max_length (int): Maximum length of sequences to include.

        Returns:
            list: A list of sequences that meet the length criteria.
        """
        # If max_length is not provided, set it to infinity
        if max_length is None:
            max_length = float('inf')

        # Filter sequences that fall within the specified length range
        filtered_sequences = [
            seq for seq in self.sequences
            if min_length <= len(seq["sequence"]) <= max_length
        ]
        # Log the number of sequences that meet the length criteria
        logging.info(f"Filtered sequences to {len(filtered_sequences)} that meet length criteria.")
        return filtered_sequences

    def subset_sequences(self, headers: list) -> list:
        """
        Subset sequences based on a list of headers.

        Args:
            headers (list): List of headers to subset sequences.

        Returns:
            list: A list of sequences that match the provided headers.
        """
        # Subset sequences by matching headers
        subset = [seq for seq in self.sequences if seq["header"] in headers]
        # Log the number of sequences that match the provided headers
        logging.info(f"Subset sequences to {len(subset)} that match provided headers.")
        return subset

    def plot_sequence_length_distribution(self, bins: int = 20, color: str = 'blue', edgecolor: str = 'black',
                                          save_path: str = None) -> None:
        """
        Plot the distribution of sequence lengths.

        Args:
            bins (int): Number of bins in the histogram.
            color (str): Color of the bars in the histogram.
            edgecolor (str): Color of the edges of the bars.
            save_path (str): Path to save the plot image file. If None, the plot is shown without saving.

        Raises:
            SequenceParserError: If matplotlib is not installed.
        """
        try:
            # Import matplotlib for plotting
            import matplotlib.pyplot as plt
            # Calculate the lengths of all sequences
            lengths = [len(seq["sequence"]) for seq in self.sequences]
            # Create a histogram of the sequence lengths with customization options
            plt.hist(lengths, bins=bins, color=color, edgecolor=edgecolor)
            # Set the title and labels of the plot
            plt.title("Sequence Length Distribution")
            plt.xlabel("Sequence Length")
            plt.ylabel("Frequency")
            # Save the plot to a file if save_path is provided
            if save_path:
                plt.savefig(save_path)
                logging.info(f"Sequence length distribution plot saved to {save_path}.")
            else:
                # Otherwise, display the plot
                plt.show()
                logging.info("Displayed sequence length distribution plot.")
        except ImportError as e:
            # Log an error if matplotlib is not installed
            logging.error("matplotlib is required to plot sequence length distribution.")
            # Raise a custom SequenceParserError if matplotlib is not installed
            raise SequenceParserError("matplotlib is required to plot sequence length distribution.") from e

    def compute_sequence_statistics(self) -> dict:
        """
        Compute statistics of the sequences.

        This method calculates and returns statistics such as the number of sequences,
        the average sequence length, and the maximum and minimum sequence lengths.

        Returns:
            dict: A dictionary containing sequence statistics.
        """
        # Calculate the lengths of all sequences
        lengths = [len(seq["sequence"]) for seq in self.sequences]
        # Calculate and return various statistics about the sequences
        stats = {
            "total_sequences": len(self.sequences),
            "average_length": sum(lengths) / len(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "min_length": min(lengths) if lengths else 0
        }
        return stats
