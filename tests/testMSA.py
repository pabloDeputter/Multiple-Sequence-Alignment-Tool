import sys
import os
import unittest

# Add src directory to sys.path
sys.path.append(os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../src')))

from msa import MSA
from utils import Utils


class TestMSA(unittest.TestCase):
    """
    Test class for the MSA module.
    """

    def setUp(self):
        """
        Set up for the test cases.
        """
        self.test_files_dir = os.path.join(
            os.path.dirname(__file__), 'data/input')
        self.expected_files_dir = os.path.join(
            os.path.dirname(__file__), 'data/expected')

        # Retrieve settings from JSON file
        settings = Utils.readSettings(os.path.abspath(os.path.join(
            os.path.dirname(__file__), '../data/settings.json')))
        self.msa = MSA(settings["match"], settings["mismatch"],
                       settings["gap_penalty"], settings["gap_gap_penalty"])

    def read_expected_output(self, file_path: str) -> str:
        """
        Read the expected output from a file.

        @param file_path: Path to the file containing the expected output.

        @return: The expected output as a string.
        """
        with open(file_path, 'r') as file:
            return file.read().strip()

    def read_input_file(self, file_path: str) -> str:
        """
        Read the input sequences from a file.

        @param file_path: Path to the file containing the input sequences.

        @return: The input sequences as a string.
        """
        with open(file_path, 'r') as file:
            return file.read().strip()


def create_test_method(file: str, alignment_type: str) -> callable:
    """
    Create a test method for the given file and alignment type.

    @param file: The name of the file to test.
    @param alignment_type: The type of alignment to perform.

    @return: The test method.
    """

    def test(self):
        """
        Test method for the given file and alignment type.

        @param self: The TestMSA instance.
        """
        test_file = os.path.join(self.test_files_dir, file)
        expected_output_file = os.path.join(self.expected_files_dir, file.replace(
            ".fasta", f".fasta_{alignment_type}_aligned.txt"))
        expected_output = self.read_expected_output(expected_output_file)
        input_sequences = self.read_input_file(test_file)

        output, _, _ = self.msa.align(test_file, True if alignment_type ==
                                                               "local" else False)

        self.assertEqual(expected_output, output)

        print(f"\nTest passed for {file} ({alignment_type} alignment):\n")
        print(f"Input sequences:\n{input_sequences}\n")
        print(f"Expected output:\n{expected_output}\n")
        print(f"Generated output:\n{output}\n")

    return test


# Dynamically create test methods for each file
for file in os.listdir(os.path.join(os.path.dirname(__file__), 'data/input')):
    if file.endswith(".fasta"):
        setattr(TestMSA, f"test_global_alignment_{
        file.split('.')[0]}", create_test_method(file, "global"))
        setattr(TestMSA, f"test_local_alignment_{
        file.split('.')[0]}", create_test_method(file, "local"))

if __name__ == '__main__':
    unittest.main()
