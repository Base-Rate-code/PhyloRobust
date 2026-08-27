"""
PhyloRobust V1
A reproducible tool for assessing the robustness of phylogenetic
conclusions across alternative distance models and tree-building methods.

Version: 1.0
"""

# ============================================================
# Imports
# ============================================================


# ============================================================
# FASTA handling
# ============================================================

def read_fasta(filename):
    """
    Read a FASTA file and return sequences as a dictionary.

    Parameters
    ----------
    filename : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Dictionary mapping sequence names to sequences.
    """

    sequences = {}
    current_name = None

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_name = line[1:].strip()

                if not current_name:
                    raise ValueError(
                        f"Duplicate sequence name found: {current_name}"
                    )

                sequences[current_name] = ""

            else:
                if current_name is None:
                    raise ValueError(
                        "Sequence data found before the first FASTA header."
                    )

                sequences[current_name] += line.upper()

    if not sequences:
        raise ValueError("No sequences found in FASTA file.")

    return sequences

def validate_alignment(sequences):
    """
    Validate that sequences form a usable multiple sequence alignment.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping sequence names to sequences.

    Returns
    -------
    bool
        True if the alignment is valid.
    """

    if len(sequences) < 2:
        raise ValueError(
            "At least two sequences are required for phylogenetic analysis."
        )

    lengths  = {len(sequence) for sequence in sequences.values()}

    if len(lengths) != 1:
        raise ValueError(
            "All sequences must have the same length."
        )

    if 0 in lengths:
        raise ValueError(
            "Sequences cannot be empty."
        )

    return True



# ============================================================
# Distance calculations
# ============================================================

def p_distance(sequence_1, sequence_2):
    """
    Calculate the p-distance between two aligned sequences.

    Gap-containing positions are excluded from the calculation.

    Parameters
    ----------
    sequence_1 : str
        First aligned sequence.
    sequence_2 : str
        Second aligned sequence.

    Returns
    -------
    float
        Proportion of comparable sites that differ.
    """

    if len(sequence_1)!= len(sequence_2):
        raise ValueError(
            "Sequences must be of same length."
        )

    differences = 0
    comparable_sites = 0

    for base_1, base_2 in zip(sequence_1, sequence_2):

        if base_1=="-" or base_2=="-":
            continue

        comparable_sites +=1

        if base_1 != base_2:
            differences +=1

    if comparable_sites==0:
        raise ValueError(
            "No comparable sites found between the two sequences."
        )

    return differences / comparable_sites



# ============================================================
# Tree construction
# ============================================================


# ============================================================
# Tree analysis
# ============================================================


# ============================================================
# Robustness analysis
# ============================================================


# ============================================================
# Output
# ============================================================


# ============================================================
# Main program
# ============================================================

def main():
    """Run PhyloRobust V1."""
    sequences = read_fasta("Ver_01/test.fasta")
    validate_alignment(sequences)

    distance = p_distance(
        sequences["Sequence_A"],
        sequences["Sequence_B"]
    )

    print("p-distance:", distance)


if __name__ == "__main__":
    main()