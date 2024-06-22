import numpy as np

from typing import List, Tuple
from tqdm import tqdm

# Own imports
from src.utils import Utils


class MSA:
    """
    Class for Multiple Sequence Alignment.
    """

    def __init__(self, match: int, mismatch: int, gap_penalty: int, gap_gap_penalty: int):
        """
        Initialize the alignment parameters.

        @param match: Score for matching characters.
        @param mismatch: Score for mismatching characters.
        @param gap_penalty: Penalty for opening a gap/indel (linear gap penalty).
        @param gap_gap_penalty: Penalty for two gaps in the same column.
        """
        self.match = match
        self.mismatch = mismatch
        self.gap_penalty = gap_penalty
        self.gap_gap_penalty = gap_gap_penalty

    def initializeScoringMatrix(self, dimensions: Tuple[int, ...], local: bool = False) -> np.ndarray:
        """
        Initialize scoring matrix for global or local alignment.

        Sequences: s1, s2, s3 with lengths 2, 2, 2 respectively (dimensions: (3, 3, 3))

        Initial 3D Matrix (all zeros):
        [
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
        ]

        After setting initial gap penalties (gap_penalty = -4):
        [
        [
            [0, -4, -8],
            [-4, -8, -12],
            [-8, -12, -16]
        ],
        [
            [-4, -8, -12],
            [-8, -12, -16],
            [-12, -16, -20]
        ],
        [
            [-8, -12, -16],
            [-12, -16, -20],
            [-16, -20, -24]
        ]
        ]

        @param dimensions: Dimensions of the scoring matrix, corresponding to length of sequence + 1.
        @param local: Whether to perform local alignment.

        @return: Initialized scoring matrix.
        """
        # Create matrix of zeros with specified dimensions
        matrix = np.zeros(dimensions, dtype=int)

        # If global alignment, values of first and rows are set to account for gap penalties.
        # For local alignment, all initial values are 0, allowing alignment from anywhere.
        if not local:
            for i in tqdm(np.ndindex(matrix.shape),
                          desc="Initializing scoring matrix",
                          total=np.prod(matrix.shape)):
                # If first cel (origin), set to 0
                if sum(i) == 0:
                    matrix[i] = 0
                # Otherwise, set to gap penalty times the sum of the indices
                else:
                    matrix[i] = sum(i) * self.gap_penalty
        return matrix

    def calculateScore(self, pred_idx: Tuple[int, ...], sequences: List[str], delta: Tuple[int, ...]) -> int:
        """
        Calculate the score for a given set of predecessor indices.

        @param pred_idx: Indices of the predecessor cell.
        @param sequences: List of sequences to be aligned.
        @param delta: Tuple indicating the movement direction for each sequence.

        @return: Column score.
        """
        # Initialize score to 0
        column_score = 0
        num_sequences = len(sequences)

        # Iterate over all pairwise combinations of sequences
        for i in range(num_sequences):
            for j in range(i + 1, num_sequences):
                # Case 1: Both sequences have moved
                if delta[i] == delta[j] == 1:
                    # If characters are a match, add match score
                    if sequences[i][pred_idx[i]] == sequences[j][pred_idx[j]]:
                        column_score += self.match
                    else:
                        # If characters are a mismatch, add mismatch score
                        column_score += self.mismatch
                # Case 2: One sequence has move, add gap penalty
                elif (delta[i] == 1 and delta[j] == 0) or (delta[i] == 0 and delta[j] == 1):
                    column_score += self.gap_penalty
                # Case 3: Both sequences have not moved, add gap-gap penalty
                elif delta[i] == delta[j] == 0:
                    column_score += self.gap_gap_penalty

        return column_score

    def getPossiblePredecessors(self, idx: Tuple[int, ...], sequences: List[str], matrix: np.ndarray) -> List[Tuple[Tuple[int, ...], int]]:
        """
        Get possible predecessor indices and their scores.

        @param idx: Current index in the scoring matrix.
        @param sequences: List of sequences to be aligned.
        @param matrix: Scoring matrix.

        @return: List of tuples containing predecessor indices and their scores.
        """
        num_sequences = len(sequences)
        possible_predecessors = []

        # Iterate over all possible movements or predecessors, in general there are 2^N - 1 possible movements
        for delta in np.ndindex((2,) * num_sequences):
            # Skip if there is no movement (delta is all zeros)
            if sum(delta) == 0:
                continue
            # Calculate indices of predecessor cell
            pred_idx = tuple(i - d for i, d in zip(idx, delta))
            # Check if all indices are valid/non-negative
            if all(i >= 0 for i in pred_idx):
                # Calculate score for the current cell based on the predecessor
                column_score = self.calculateScore(
                    pred_idx, sequences, delta)
                # Add predecessor indices and score to possible predecessors
                possible_predecessors.append(
                    (pred_idx, matrix[pred_idx] + column_score))

        return possible_predecessors

    def fillScoringMatrix(self, matrix: np.ndarray, sequences: List[str], local: bool = False) -> np.ndarray:
        """
        Fill in the scoring matrix, calculate alignment scores for either global or local alignment.

        @param matrix: Scoring matrix.
        @param sequences: List of sequences to be aligned.
        @param local: Whether to perform local alignment.

        @return: Scoring matrix with calculated scores.
        """
        dims = matrix.shape

        for idx in tqdm(np.ndindex(dims), desc="Filling scoring matrix", total=np.prod(dims)):
            # Skip origin cell (0, 0, ..., 0)
            if sum(idx) == 0:
                continue
            # Local alignment can start from anywhere hence 0
            possible_scores = [0] if local else []
            # Get possible predecessors and their scores
            possible_predecessors = self.getPossiblePredecessors(
                idx, sequences, matrix)
            # Add scores of possible predecessors to possible scores
            possible_scores.extend(score for _, score in possible_predecessors)
            # Set the current cell to the maximum score from the possible scores
            matrix[idx] = max(possible_scores)

        return matrix

    def traceback(self, matrix: np.ndarray, sequences: List[str], local: bool = False) -> Tuple[List[str], int]:
        """
        Perform traceback to reconstruct the aligned sequences from the scoring matrix.

        @param matrix: Scoring matrix.
        @param sequences: List of sequences to be aligned.
        @param local: Whether to perform local alignment.

        @return: A tuple containing the aligned sequences and the maximum score.
        """
        # Initialize list of aligned sequences
        aligned_sequences = ['' for _ in sequences]

        # If local, find the maximum score and its index
        if local:
            # Find index of maximum value in the matrix and convert to tuple
            idx = np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)
            max_score = matrix[idx]
        else:
            # If global, set index to the last cell in the matrix
            idx = tuple(dim - 1 for dim in matrix.shape)
            max_score = matrix[idx]

        # Perform traceback while: (1) for local alignment, the score is greater than 0, or (2) for global
        # alignment, any index is greater than 0. So we stop when either the score is 0 or we reach the origin cell.
        while (local and matrix[idx] > 0) or (not local and any(i > 0 for i in idx)):
            # Get possible predecessors indices and their scores
            possible_predecessors = self.getPossiblePredecessors(
                idx, sequences, matrix)
            # Find the predecessor with the maximum score
            best_pred_idx, _ = max(
                possible_predecessors, key=lambda x: x[1])
            # Calculate delta between current index and predecessor index
            delta = tuple(i - j for i, j in zip(idx, best_pred_idx))
            idx = best_pred_idx
            for i, d in enumerate(delta):
                # If there is movement in the sequence, add the character to the aligned sequence
                if d:
                    aligned_sequences[i] = sequences[i][idx[i]] + \
                        aligned_sequences[i]
                else:
                    # If there is no movement, add a gap
                    aligned_sequences[i] = '.' + aligned_sequences[i]

        return aligned_sequences, max_score

    def align(self, fasta_file: str, local: bool = False) -> Tuple[str, List[str], int]:
        """
        Perform multiple sequence alignment.

        @param fasta_file: Path to the FASTA file containing the sequences.
        @param local: Whether to perform local alignment.

        @return: A tuple containing the aligned sequences, the alignment score, and the final score.
        """
        sequences, ids = Utils.parseFasta(fasta_file)
        # Calculate dimensions of scoring matrix, add 1 for initial row and
        # column for gaps at beginning of sequences. For example:
        """
        >s1
        ACTG
        >s2
        TG
        gives dimensions (4, 3) and scoring matrix:
        -  A  C  T  G
        -  0 -4 -8-12-16
        T -4 -2 -6 -3 -7
        G -8 -6 -4 -7 2
        resulting in alignment with score 2:
        s1: ACTG
        s2: ..TG
        """
        dimensions = tuple(len(seq) + 1 for seq in sequences)

        # Initialize scoring matrix
        score_matrix = self.initializeScoringMatrix(dimensions, local=local)
        # Fill scoring matrix
        score_matrix = self.fillScoringMatrix(
            score_matrix, sequences, local=local)
        print("Performing traceback...\n")
        # Perform traceback
        aligned_sequences, final_score = self.traceback(
            score_matrix, sequences, local=local)
        
        output = []
        for seq_id, aligned_seq in zip(ids, aligned_sequences):
            output.append(f"{seq_id}: {aligned_seq}")
        output.append(f"Alignment score: {final_score}")
        output = "\n".join(output)

        return output, aligned_sequences, final_score


if __name__ == "__main__":
    # Initialize MSA
    settings_file = "../data/settings.json"
    settings = Utils.readSettings(settings_file)
    fasta_file = "../tests/data/input/test6.fasta"
    msa = MSA(settings["match"], settings["mismatch"], settings["gap_penalty"], settings["gap_gap_penalty"])

    output, aligned_sequences, final_score = msa.align(fasta_file, local=True)
    print(output)
