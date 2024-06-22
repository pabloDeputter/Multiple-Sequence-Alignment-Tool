import argparse
import json

from src.msa import MSA
from src.utils import Utils

# Default settings: match = 5, mismatch = -2, gap_penalty = -4, gap_gap_penalty = 0
DEFAULT_SETTINGS_PATH = "data/settings.json"


def main():
    """
    Main function for the CLI. Can be run:
        - With no settings provided (uses default settings):
            > python main.py -i data/input/input.fasta -o data/output/output.txt -a global
        - With settings provided (in correct format, otherwise default settings are used):
            > python main.py -i data/input/input.fasta -o data/output/output.txt -s data/settings.json -a global
    """
    parser = argparse.ArgumentParser(
        description="Multiple Sequence Alignment Tool")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to the input FASTA file")
    parser.add_argument("-o", "--output", required=False,
                        help="Path to the output file")
    parser.add_argument("-s", "--settings",
                        help="Path to the settings JSON file")
    parser.add_argument("-a", "--alignment", choices=[
                        "local", "global"], default="global", help="Alignment type (default: global)")

    args = parser.parse_args()

    # Read settings from JSON file if provided, otherwise use default settings
    if args.settings:
        try:
            settings = Utils.readSettings(args.settings)
        except (FileNotFoundError, json.JSONDecodeError, ValueError) as e:
            print(f"Error loading settings file: {e}")
            print("Using default settings instead")
            settings = Utils.readSettings(DEFAULT_SETTINGS_PATH)
    else:
        settings = Utils.readSettings(DEFAULT_SETTINGS_PATH)

    # Initialize MSA with settings
    msa = MSA(settings["match"], settings["mismatch"],
              settings["gap_penalty"], settings["gap_gap_penalty"])

    # Perform alignment based on the specified alignment type
    if args.alignment == "local":
        output, _, _ = msa.align(
            args.input, local=True)
    else:
        output, _, _ = msa.align(
            args.input, local=False)

    if args.output:
        # Write output to file
        with open(args.output, "w") as f:
            f.write(output)

    print(output)


if __name__ == "__main__":
    # Initialize MSA
    # settings_file = "data/settings.json"
    # settings = Utils.readSettings(settings_file)
    # fasta_file = "tests/data/input/test4.fasta"
    # msa = MSA(settings["match"], settings["mismatch"], settings["gap_penalty"], settings["gap_gap_penalty"])
    #
    # output, aligned_sequences, final_score = msa.align(fasta_file, local=True)
    # print(output)

    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align.Applications import ClustalwCommandline

    # Align the sequences using ClustalW
    clustalw_cline = ClustalwCommandline("clustalw2", infile="cs_assignment.fasta")
    stdout, stderr = clustalw_cline()
    alignment = AlignIO.read("cs_assignment.aln", "clustal")
    print(alignment)
    summary_align = AlignInfo.SummaryInfo(alignment)
    pssm = summary_align.pos_specific_score_matrix()
    print(pssm)
    conserved_region = ""
    threshold = 0.7  # Example threshold for conservation
    for position in pssm:
        max_score = max(position.values())
        if max_score > threshold:
            conserved_region += max(position, key=position.get)
        else:
            conserved_region += "-"
    print("Conserved Region:", conserved_region)


    # CLI
    # main()

