# bioseqlab/parsers/sequence_parser.py
from abc import ABC, abstractmethod
import logging
import os


# Define a custom exception for handling specific sequence parser errors
class SequenceParserError(Exception):
    """Custom exception for errors in the SequenceParser class."""
    pass


# Define the abstract base class for sequence parsing
class SequenceParser(ABC):
    """
    Abstract base class for parsing sequence files.

    This class provides a common interface and shared functionality for all sequence parsers.
    Subclasses must implement the abstract methods defined here.
    """

    def __init__(self, filename: str, log_level: int = logging.INFO):
        """
        Initialize the SequenceParser with the input filename.

        The constructor initializes the class with a filename and sets up logging.
        It also immediately validates the provided file to ensure it exists and is accessible.

        Args:
            filename (str): The path to the sequence file.
            log_level (int): The logging level (default is INFO).

        Raises:
            SequenceParserError: If the filename is not provided or the file is invalid.
        """
        # Set up logging configuration to handle different levels of logging output
        logging.basicConfig(
            level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

        # Store the filename provided during initialization in an instance variable
        self.filename: str = filename

        # Initialize an empty list to hold the parsed sequences
        self.sequences: list = []

        # Call the method to validate the input file upon initialization
        self._validate_file()

    def _validate_file(self) -> None:
        """
        Validate that the input file exists, is accessible, and has the correct format.

        This method checks whether the file exists, is accessible, and can be opened.
        If any of these checks fail, it raises a custom SequenceParserError with a descriptive message.

        Raises:
            SequenceParserError: If the file does not exist, is not accessible, or is not a recognized format.
        """
        # Check if the filename was provided
        if not self.filename:
            logging.error("No filename provided.")
            raise SequenceParserError("No filename provided.")

        # Check if the file exists and is accessible in the filesystem
        if not os.path.isfile(self.filename):
            logging.error(f"File not found: {self.filename}")
            raise SequenceParserError(f"File not found: {self.filename}")

        # Attempt to open the file to ensure it is readable
        try:
            with open(self.filename, 'r') as file:
                logging.info(f"Successfully opened file: {self.filename}")
        except IOError as e:
            # Handle IO errors and raise a custom SequenceParserError
            logging.error(
                f"IO error while accessing the file: {self.filename}")
            raise SequenceParserError(
                f"IO error while accessing the file: {self.filename}") from e

    def _validate_file_type(self, allowed_extensions: list) -> None:
        """
        Validate the file type based on its extension.

        This method ensures that the file being processed has a valid extension.
        It checks the file extension against a list of allowed extensions and raises an error if the extension is not valid.

        Args:
            allowed_extensions (list): List of allowed file extensions (e.g., ['.fasta', '.fa']).

        Raises:
            SequenceParserError: If the file extension is not allowed.
        """
        # Extract the file extension from the filename
        file_extension = os.path.splitext(self.filename)[1]

        # Check if the file extension is in the list of allowed extensions
        if file_extension not in allowed_extensions:
            logging.error(f"File extension {file_extension} is not supported.")
            raise SequenceParserError(
                f"Unsupported file extension: {file_extension}")

    def _validate_sequences(self) -> None:
        """
        Validate the parsed sequences to ensure they meet expected criteria.

        This method can be extended by subclasses to include specific validation
        rules for different file formats. Additionally it checks whether any sequences were successfully parsed.
        If no sequences are found, it logs a warning and raises a custom SequenceParserError.

        Raises:
            SequenceParserError: If no sequences are found or sequences are improperly formatted.
        """
        # Check if any sequences have been parsed
        if not self.sequences:
            logging.warning(
                "No sequences found. Please ensure the file is correctly formatted and parsed.")
            raise SequenceParserError(
                "No sequences found. Please ensure the file is correctly formatted and parsed.")
        else:
            # Log the number of validated sequences
            logging.info(f"Validated {len(self.sequences)} sequences.")

    @property
    @abstractmethod
    def supported_extensions(self) -> list:
        """
        Abstract property that returns the list of supported file extensions for the parser.

        This property enforces that each subclass specifies which file extensions it supports.
        It must be implemented by each subclass to return a list of valid extensions (e.g., ['.fasta', '.fa'] for FASTA files).

        Returns:
            list: List of supported file extensions.

        Raises:
            NotImplementedError: If the subclass does not implement this property.
        """
        pass

    @abstractmethod
    def read_sequences(self) -> None:
        """
        Abstract method to read sequences from the file.

        Subclasses must implement this method to handle the specific logic for reading sequences from the file format they support.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def get_headers(self) -> list:
        """
        Abstract method to retrieve the headers of the sequences.

        This method should be implemented by subclasses to return a list of headers for the sequences in the file.
        The headers typically contain identifier information for each sequence.

        Returns:
            list: A list of sequence headers.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def get_sequence_by_header(self, header: str) -> str:
        """
        Abstract method to retrieve a sequence by its header.

        This method should be implemented by subclasses to return the sequence corresponding to a given header.
        This allows users to retrieve specific sequences by referencing their identifiers.

        Args:
            header (str): The header of the sequence.

        Returns:
            str: The sequence corresponding to the header.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def plot_sequence_length_distribution(self) -> None:
        """
        Abstract method to plot the distribution of sequence lengths.

        This method should be implemented by subclasses to create a plot that shows the distribution of sequence lengths.
        It is useful for visualizing the variability in sequence sizes within the file.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def compute_sequence_statistics(self) -> dict:
        """
        Abstract method to compute statistics of sequences.

        This method should be implemented by subclasses to calculate and return statistics about the sequences, such as average length, GC content, etc.
        The statistics can provide insights into the characteristics of the sequences in the file.

        Returns:
            dict: A dictionary containing sequence statistics.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        This method allows the class to be used with a 'with' statement, ensuring proper resource management.

        Returns:
            SequenceParser: The object itself, allowing it to be used in a 'with' statement.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit the runtime context related to this object.

        This method handles any cleanup actions needed upon exiting the 'with' statement.
        It can be used to perform cleanup tasks, such as closing files or releasing resources.

        Args:
            exc_type: The exception type, if an exception occurred.
            exc_val: The exception value, if an exception occurred.
            exc_tb: The traceback, if an exception occurred.
        """
        pass