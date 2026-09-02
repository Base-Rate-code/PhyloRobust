"""
PhyloRobust V1
A reproducible tool for assessing the robustness of phylogenetic
conclusions across alternative distance models and tree-building methods.

Version: 1.0
"""

# ============================================================
# Imports
# ============================================================

import math

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
                        "FASTA header cannot be empty."
                        )

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

    def __init__(self, name=None, branch_length=0.0, height = 0.0):
        self.name =name
        self.branch_length = branch_length
        self.height =  height
        self.children =[]
        self.parent = None

    def add_child(self, child):
        """
        Add a child node and assign this node as its parent.
        """

        self.children.append(child)
        child.parent = self

def find_closest_clusters(distance_matrix, clusters):
    """
    Find the pair of clusters with the smallest distance.

    Parameters
    ----------
    distance_matrix : dict
        Current distance matrix.

    clusters : dict
        Current active clusters.

        Returns
    -------
    tuple
        Names of the two closest clusters.
    """

    cluster_names = list(clusters.keys())

    best_pair = None
    best_distance = float("inf")

    for i in range(len(cluster_names)):
        for j in range(i+1, len(cluster_names)):

            cluster_1 = cluster_names[i]
            cluster_2 = cluster_names[j]

            distance = distance_matrix[cluster_1][cluster_2]

            if distance < best_distance:
                best_distance= distance
                best_pair = (cluster_1, cluster_2)

    return best_pair

def merge_clusters(clusters, cluster_1, cluster_2, distance):
    """
    Merge two clusters into a new internal node.

    Parameters
    ----------
    clusters : dict
        Current active clusters.

    cluster_1 : str
        Name of the first cluster.

    cluster_2 : str
        Name of the second cluster.

    Returns
    -------
    tuple
        New cluster name and new Node object.
    """

    node_1 = clusters[cluster_1]
    node_2 = clusters[cluster_2]

    new_name = f"({cluster_1}, {cluster_2})"
    new_height= distance/ 2

    # Calculate branch lengths from the new node
    # to each child.
    node_1.branch_length = new_height - node_1.height
    node_2.branch_length = new_height - node_2.height

    new_node = Node(
        new_name,
        height=  new_height
        )

    new_node.add_child(node_1)
    new_node.add_child(node_2)

    return new_name, new_node

def update_distance_matrix(distance_matrix, cluster_1, cluster_2, new_cluster):
    """
    Update the distance matrix after merging two clusters.

    UPGMA calculates the distance between the new cluster and
    each remaining cluster using the arithmetic mean.

    Parameters
    ----------
    distance_matrix : dict
        Current distance matrix.

    cluster_1 : str
        First cluster being merged.

    cluster_2 : str
        Second cluster being merged.

    new_cluster : str
        Name of the newly created cluster.

    Returns
    -------
    dict
        Updated distance matrix.
    """

    remaining_clusters= [
        name for name in distance_matrix
        if name not in (cluster_1, cluster_2)
    ]
    new_matrix = {}

    # Create entries for the remaining clusters
    for name in remaining_clusters:
        new_matrix[name] = {}

    # Create entry for the new cluster
    new_matrix[new_cluster] = {}

    # Calculate distances between the new cluster
    # and every remaining cluster.
    for name in remaining_clusters:

        distance_1 = distance_matrix[cluster_1][name]
        distance_2 = distance_matrix[cluster_2][name]

        new_distance = (distance_1 + distance_2)/2

        new_matrix[name][new_cluster] = new_distance
        new_matrix[new_cluster][name] = new_distance

    # Distance of a cluster to itself is zero.
    new_matrix[new_cluster][new_cluster] = 0.0

    # Preserve zero self-distances for remaining clusters.
    for name in remaining_clusters:
        new_matrix[name][name] = 0.0

    return new_matrix

def upgma (distance_matrix, clusters):
    """
    Construct a phylogenetic tree using UPGMA.

    Parameters
    ----------
    distance_matrix : dict
        Initial pairwise distance matrix.

    clusters : dict
        Dictionary mapping sequence names to Node objects.

    Returns
    -------
    Node
        Root node of the UPGMA tree.
    """
    while len(clusters) > 1:

        # Find the closest pair of clusters.
        cluster_1, cluster_2 = find_closest_clusters(
            distance_matrix,
            clusters
        )

        distance = distance_matrix[cluster_1][cluster_2]

        # Merge the closest clusters.
        new_name, new_node = merge_clusters(
            clusters,
            cluster_1,
            cluster_2,
            distance
        )


        # Update the distance matrix.
        distance_matrix = update_distance_matrix(
            distance_matrix,
            cluster_1,
            cluster_2,
            new_name
        )

        # Remove the old clusters.
        del clusters[cluster_1]
        del clusters[cluster_2]

        # Add the new cluster.
        clusters[new_name] = new_node

    # The final remaining cluster is the root.
    root = next(iter(clusters.values()))

    return root


def calculate_row_totals(distance_matrix):
    """
    Calculate the total distance from each cluster
    to all other clusters.
    """

    row_totals = {}

    for name, distances in distance_matrix.items():

        total = sum(
            distance
            for other_name, distance in distances.items()
            if other_name != name
        )

        row_totals[name] = total

    return row_totals


def calculate_q_matrix(distance_matrix):
    """
    Calculate the Neighbor Joining Q-matrix.

    Q(i,j) = (n-2)d(i,j) - r(i) - r(j)
    """

    names = list(distance_matrix.keys())
    n = len(names)

    row_totals = calculate_row_totals(distance_matrix)

    q_matrix = {}

    for name_1 in names:
        q_matrix[name_1] = {}

        for name_2 in names:

            if name_1 == name_2:
                q_matrix[name_1][name_2] = 0.0

            else:
                distance = distance_matrix[name_1][name_2]

                q_value = (
                    (n - 2) * distance
                    - row_totals[name_1]
                    - row_totals[name_2]
                )

                q_matrix[name_1][name_2] = q_value

    return q_matrix


def find_minimum_q_pair(q_matrix):
    """
    Find the pair of clusters with the smallest Q-value.
    """

    names = list(q_matrix.keys())

    best_pair = None
    best_value = float("inf")

    for i in range(len(names)):

        for j in range(i + 1, len(names)):

            name_1 = names[i]
            name_2 = names[j]

            q_value = q_matrix[name_1][name_2]

            if q_value < best_value:

                best_value = q_value
                best_pair = (name_1, name_2)

    return best_pair


def calculate_nj_branch_lengths(
    distance_matrix,
    row_totals,
    cluster_1,
    cluster_2
):
    """
    Calculate the Neighbor Joining branch lengths
    from the new internal node to the two selected clusters.
    """

    n = len(distance_matrix)

    distance = distance_matrix[cluster_1][cluster_2]

    branch_1 = (
        0.5 * distance
        + (row_totals[cluster_1] - row_totals[cluster_2])
        / (2 * (n - 2))
    )

    branch_2 = distance - branch_1

    return branch_1, branch_2


def create_nj_node(
    clusters,
    cluster_1,
    cluster_2,
    branch_1,
    branch_2
):
    """
    Create an internal Neighbor Joining node.
    """

    node_1 = clusters[cluster_1]
    node_2 = clusters[cluster_2]

    node_1.branch_length = branch_1
    node_2.branch_length = branch_2

    new_name = f"({cluster_1}, {cluster_2})"

    new_node = Node(new_name)

    new_node.add_child(node_1)
    new_node.add_child(node_2)

    return new_name, new_node


def update_nj_distance_matrix(
    distance_matrix,
    cluster_1,
    cluster_2,
    new_cluster
):
    """
    Reduce the Neighbor Joining distance matrix
    after joining two clusters.

    d(u,k) = [d(i,k) + d(j,k) - d(i,j)] / 2
    """

    remaining_clusters = [
        name
        for name in distance_matrix
        if name not in (cluster_1, cluster_2)
    ]

    new_matrix = {}

    for name in remaining_clusters:
        new_matrix[name] = {}

    new_matrix[new_cluster] = {}

    for name in remaining_clusters:

        distance_1 = distance_matrix[cluster_1][name]
        distance_2 = distance_matrix[cluster_2][name]
        distance_12 = distance_matrix[cluster_1][cluster_2]

        new_distance = (
            distance_1
            + distance_2
            - distance_12
        ) / 2

        new_matrix[name][new_cluster] = new_distance
        new_matrix[new_cluster][name] = new_distance

    for name in remaining_clusters:
        new_matrix[name][name] = 0.0

    new_matrix[new_cluster][new_cluster] = 0.0

    return new_matrix


def nj(distance_matrix, clusters):
    """
    Construct a phylogenetic tree using Neighbor Joining.
    """

    while len(clusters) > 2:

        # Step 1: calculate row totals
        row_totals = calculate_row_totals(distance_matrix)

        # Step 2: calculate Q-matrix
        q_matrix = calculate_q_matrix(distance_matrix)

        # Step 3: find the pair with minimum Q
        cluster_1, cluster_2 = find_minimum_q_pair(q_matrix)

        # Step 4: calculate branch lengths
        branch_1, branch_2 = calculate_nj_branch_lengths(
            distance_matrix,
            row_totals,
            cluster_1,
            cluster_2
        )
        

        # Step 5: create internal node
        new_name, new_node = create_nj_node(
            clusters,
            cluster_1,
            cluster_2,
            branch_1,
            branch_2
        )

        # Step 6: reduce distance matrix
        distance_matrix = update_nj_distance_matrix(
            distance_matrix,
            cluster_1,
            cluster_2,
            new_name
        )

        # Step 7: remove old clusters
        del clusters[cluster_1]
        del clusters[cluster_2]

        # Step 8: add new cluster
        clusters[new_name] = new_node

    # Final two clusters
    remaining = list(clusters.keys())

    cluster_1 = remaining[0]
    cluster_2 = remaining[1]

    distance = distance_matrix[cluster_1][cluster_2]

    node_1 = clusters[cluster_1]
    node_2 = clusters[cluster_2]

    node_1.branch_length = distance / 2
    node_2.branch_length = distance / 2

    root = Node(
        f"({cluster_1}, {cluster_2})"
    )

    root.add_child(node_1)
    root.add_child(node_2)

    return root

def run_pipeline(sequences, distance_function, tree_method):
    """
    Run one phylogenetic analysis pipeline.

    Parameters
    ----------
    sequences : dict
        Aligned sequences.

    distance_function : function
        Distance model to use.

    tree_method : function
        Tree-building method to use.

    Returns
    -------
    Node
        Root node of the resulting phylogenetic tree.
    """

    # Calculate pairwise distances.
    distance_matrix = calculate_distance_matrix(
        sequences,
        distance_function
    )

    # Create a fresh set of leaf nodes.
    clusters = {
        name: Node(name)
        for name in sequences
    }

    # Build the tree.
    tree = tree_method(
        distance_matrix,
        clusters
    )

    return tree


def get_leaf_names(node):
    """
    Return the names of all leaf sequences below a node.
    """

    # If the node has no children, it is a leaf.
    if not node.children:
        return [node.name]

    leaf_names = []

    for child in node.children:
        leaf_names.extend(
            get_leaf_names(child)
        )

    return leaf_names

def get_clades(node):
    """
    Return all non-trivial clades in a tree.

    A clade is represented as a frozenset of leaf names.
    """

    clades = []

    if not node.children:
        return clades

    leaves = frozenset(get_leaf_names(node))

    if len(leaves) > 1:
        clades.append(leaves)

    for child in node.children:
        clades.extend(get_clades(child))

    return clades

def clade_supported(tree, target_group):
    """
    Check whether a target group of sequences
    occurs as a clade in the tree.
    """

    target_group = frozenset(target_group)

    clades = get_clades(tree)

    for clade in clades:

        if clade == target_group:
            return True

    return False

def calculate_support(tree_results, target_group):
    """
    Calculate how many phylogenetic trees support
    a target biological grouping.
    """

    supported = 0
    total = len(tree_results)

    for tree in tree_results:
        if clade_supported(tree, target_group):
            supported += 1

    support_fraction = supported / total

    return supported, total, support_fraction







# ============================================================
# Tree analysis
# ============================================================

def print_tree(node, level=0):
    """
    Print the tree structure for inspection.
    """

    print(
        "  " * level,
        node.name,
        "height =", node.height,
        "branch =", node.branch_length
    )

    for child in node.children:
        print_tree(child, level + 1)

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

    # --------------------------------------------------
    # 1. Read input sequences
    # --------------------------------------------------

    sequences = read_fasta("Ver_01/test.fasta")

    # --------------------------------------------------
    # 2. Validate the alignment
    # --------------------------------------------------

    validate_alignment(sequences)

    # --------------------------------------------------
    # 3. Display the p-distance matrix
    # --------------------------------------------------

    matrix = calculate_distance_matrix(
        sequences,
        p_distance
    )

    print("Distance matrix:")

    for name, distances in matrix.items():
        print(name, distances)

    # --------------------------------------------------
    # 4. Run all four phylogenetic pipelines
    # --------------------------------------------------

    p_distance_upgma = run_pipeline(
        sequences,
        p_distance,
        upgma
    )

    p_distance_nj = run_pipeline(
        sequences,
        p_distance,
        nj
    )

    jc_upgma = run_pipeline(
        sequences,
        jukes_cantor_distance,
        upgma
    )

    jc_nj = run_pipeline(
        sequences,
        jukes_cantor_distance,
        nj
    )

        # --------------------------------------------------
    # 5. Define the biological claim being tested
    # --------------------------------------------------

    target_group = {
        "Sequence_A",
        "Sequence_B"
    }

    # --------------------------------------------------
    # 6. Collect all pipeline results
    # --------------------------------------------------

    tree_results = [
        p_distance_upgma,
        p_distance_nj,
        jc_upgma,
        jc_nj
    ]

    # --------------------------------------------------
    # 7. Calculate robustness
    # --------------------------------------------------

    supported, total, fraction = calculate_support(
        tree_results,
        target_group
    )

    print("\nRobustness:")
    print(f"Target group: {target_group}")
    print(f"Supported: {supported}/{total}")
    print(f"Robustness: {fraction:.2%}")

if __name__ == "__main__":
    main()
