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

def calculate_distance_matrix(sequences, distance_function):
    """
    Calculate a pairwise distance matrix.

    Parameters
    ----------
    sequences : dict
        Dictionary mapping sequence names to sequences.

    distance_function : function
        Function used to calculate the distance between two sequences.

    Returns
    -------
    dict
        Nested dictionary containing pairwise distances.
    """

    names = list(sequences.keys())
    matrix={}

    for name_1 in names:
        matrix[name_1]= {}

        for name_2 in names:

            if name_1 == name_2:
                matrix[name_1][name_2] = 0.0

            else:
                matrix[name_1][name_2] = distance_function(
                    sequences[name_1],
                    sequences[name_2]
                )

    return matrix

import math

def jukes_cantor_distance(sequence_1, sequence_2):
     """
    Calculate the Jukes-Cantor evolutionary distance
    between two aligned nucleotide sequences.

    Gap-containing positions are excluded.

    Parameters
    ----------
    sequence_1 : str
        First aligned sequence.

    sequence_2 : str
        Second aligned sequence.

    Returns
    -------
    float
        Jukes-Cantor distance.
    """

     p = p_distance(sequence_1, sequence_2)

     if p >= 0.75:
        raise ValueError(
             "Jukes Cantor distanc is undefined when"
             "p-distance is >= 0.75."
         )

     return -0.75 *math.log(1-(4/3)*p)

     




# ============================================================
# Tree construction
# ============================================================

class Node:
    """
    Represent a node in a phylogenetic tree.
    """

    def __init__(self, name=None, branch_length=0.0):
        self.name =name
        self.branch_length = branch_length
        self.children =[]
        self.parent = None

    def add_child(self, child):
        """
        Add a child node and assign this node as its parent.
        """

        self.children.append(child)
        child.parent = self




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


    node_a = Node("Sequence_A")
    node_b = Node("Sequence_B")

    parent = Node("Internal")

    parent.add_child(node_a)
    parent.add_child(node_b)

    print("Parent:", parent.name)
    print("Children:")

    for child in parent.children:
        print(child.name)
        
    print("A's parent:", node_a.parent.name)
    print("B's parent:", node_b.parent.name)



if __name__ == "__main__":
    main()