import random
import json
import numpy as np

from typing import List, Tuple
from Bio import SeqIO


class Utils:
    """
    Utility class for common functions.
    """
    @staticmethod
    def readSettings(file: str) -> dict:
        """
        Read settings from a JSON file, and validate the format.
        Must contain the following: match, mismatch, gap_penalty, gap_gap_penalty.

        @param file: Path to the JSON file.

        @return: A dictionary containing the settings.
        """
        with open(file, "r") as f:
                settings = json.load(f)
        # Validate settings
        required_keys = ["match", "mismatch", "gap_penalty", "gap_gap_penalty"]
        if not all(key in settings for key in required_keys):
            raise ValueError("Invalid settings file format")
        return settings

    @staticmethod
    def parseFasta(file: str) -> Tuple[List[str], List[str]]:
        """
        Parse a FASTA file.

        @param file: Path to the FASTA file.

        @return: A tuple containing a list of sequences and a list of IDs.
        """
        sequences = []
        ids = []
        for record in SeqIO.parse(file, "fasta"):
            # Append the sequence and ID to their respective lists
            sequences.append(str(record.seq))
            ids.append(record.id)
        return sequences, ids

    @staticmethod
    def writeFasta(file: str, sequences: List[str], ids: List[str]) -> None:
        """
        Write a list of sequences to a FASTA file.

        @param file: Path to the output FASTA file.
        @param sequences: List of sequences.
        @param ids: List of IDs.
        """
        with open(file, "w") as f:
            for i in range(len(sequences)):
                f.write(f">{ids[i]}\n")
                f.write(f"{sequences[i]}\n")

    @staticmethod
    def generateRandomSequence(n: int) -> str:
        """
        Generate a random DNA sequence of length n.

        @param n: Length of the sequence.

        @return: A random DNA sequence of length n.
        """
        bases = ['A', 'C', 'G', 'T']
        return ''.join(random.choices(bases, k=n))

    @staticmethod
    def generateTestFile(file: str, num_sequences: int = 2, length_mean: int = 100, length_std: int = 20) -> None:
        """
        Generate a test FASTA file with random sequences.

        @param file: Path to the output FASTA file.
        @param num_sequences: Number of sequences to generate.
        @param length_mean: Mean length of the sequences.
        @param length_std: Standard deviation of the sequence lengths.
        """
        # Sample sequences from a normal distribution with specified mean and std
        sequences = [Utils.generateRandomSequence(
            int(np.random.normal(length_mean, length_std))) for _ in range(num_sequences)]
        ids = [f"seq{i}" for i in range(num_sequences)]
        Utils.writeFasta(file, sequences, ids)
