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
import os
import csv
import argparse

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

                if current_name in sequences:
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
def validate_target_group(sequences, target_group):
    """
    Validate that all taxa in the target group are present
    in the input alignment.

    Parameters
    ----------
    sequences : dict
        Input sequences.

    target_group : set
        Taxa defining the target clade.

    Raises
    ------
    ValueError
        If the target group is empty or contains taxa
        not present in the alignment.
    """

    if not target_group:
        raise ValueError(
            "Target group cannot be empty."
        )

    if len(target_group) < 2:
        raise ValueError(
        "Target group must contain at least two taxa."
        )

    missing = target_group - set(sequences.keys())

    if missing:
        raise ValueError(
            f"Target taxa not found in alignment: "
            f"{', '.join(sorted(missing))}"
        )

def parse_arguments():
    """
    Parse command-line arguments for PhyloRobust.
    """

    parser = argparse.ArgumentParser(
        description=(
            "PhyloRobust V1: assess the robustness of "
            "a target phylogenetic clade across "
            "multiple analytical pipelines."
        )
    )

    parser.add_argument(
        "fasta",
        help="Path to the aligned nucleotide FASTA file."
    )

    parser.add_argument(
        "--target",
        nargs="+",
        required=True,
        help=(
            "Taxa defining the target clade. "
            "Provide two or more sequence names."
        )
    )

    return parser.parse_args()
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

def tree_to_newick(node):
    """
    Convert a phylogenetic tree into Newick format.

    Parameters
    ----------
    node : Node
        Root node of the tree.

    Returns
    -------
    str
        Tree represented in Newick format.
    """

    if not node.children:
        return f"{node.name}:{node.branch_length}"

    children = [
        tree_to_newick(child)
        for child in node.children
    ]

    return f"({','.join(children)}):{node.branch_length}"

def export_newick(tree):
    """
    Convert a tree into a complete Newick string.
    """

    newick = tree_to_newick(tree)

    # The root branch length is normally zero
    # and is not required in the final Newick representation.
    if newick.endswith(":0.0"):
        newick = newick[:-4]

    return newick + ";"

def save_newick(tree, filename):
    """
    Save a phylogenetic tree in Newick format.

    Parameters
    ----------
    tree : Node
        Root node of the phylogenetic tree.

    filename : str
        Output filename.
    """

    newick = export_newick(tree)

    with open(filename, "w") as file:
        file.write(newick + "\n")

def find_negative_branches(node, warnings=None):
    """
    Find negative branch lengths in a phylogenetic tree.

    Parameters
    ----------
    node : Node
        Root node of the tree.

    warnings : list, optional
        List used to store warning messages.

    Returns
    -------
    list
        List of warning messages for negative branch lengths.
    """

    if warnings is None:
        warnings = []

    if node.branch_length < 0:
        warnings.append(
            f"Negative branch length detected for "
            f"{node.name}: {node.branch_length}"
        )

    for child in node.children:
        find_negative_branches(child, warnings)

    return warnings

# ============================================================
# Robustness analysis
# ============================================================

PIPELINES = [
    {
        "name": "p-distance + UPGMA",
        "filename": "p_distance_UPGMA",
        "distance_model": "p-distance",
        "distance_function": p_distance,
        "tree_method": "UPGMA",
        "tree_function": upgma
    },
    {
        "name": "p-distance + NJ",
        "filename": "p_distance_NJ",
        "distance_model": "p-distance",
        "distance_function": p_distance,
        "tree_method": "NJ",
        "tree_function": nj
    },
    {
        "name": "Jukes-Cantor + UPGMA",
        "filename": "jukes_cantor_UPGMA",
        "distance_model": "Jukes-Cantor",
        "distance_function": jukes_cantor_distance,
        "tree_method": "UPGMA",
        "tree_function": upgma
    },
    {
        "name": "Jukes-Cantor + NJ",
        "filename": "jukes_cantor_NJ",
        "distance_model": "Jukes-Cantor",
        "distance_function": jukes_cantor_distance,
        "tree_method": "NJ",
        "tree_function": nj
    }
]

def analyse_pipeline(sequences, pipeline, target_group):

    tree = run_pipeline(
        sequences,
        pipeline["distance_function"],
        pipeline["tree_function"]
    )

    supported = clade_supported(
        tree,
        target_group
    )

    warnings = find_negative_branches(tree)

    if supported:
        conclusion = "SUPPORTED"
    else:
        conclusion = "REJECTED"

    result = {
        "pipeline": pipeline["name"],
        "filename": pipeline["filename"],
        "distance_model": pipeline["distance_model"],
        "tree_method": pipeline["tree_method"],
        "target_group": target_group,
        "target_supported": supported,
        "conclusion": conclusion,
        "warnings": warnings,
        "tree": tree
    }
    return result

def calculate_robustness(pipeline_results):
    """
    Calculate overall robustness across all tested pipelines.

    Parameters
    ----------
    pipeline_results : list
        Results generated by analyse_pipeline().

    Returns
    -------
    dict
        Summary of supported, rejected, total pipelines,
        robustness, and warnings.
    """

    total = len(pipeline_results)

    if total == 0:
        return {
            "supported": 0,
            "rejected": 0,
            "total": 0,
            "robustness": 0.0,
            "warning_pipelines": 0
        }

    supported = sum(
        result["target_supported"]
        for result in pipeline_results
    )

    rejected = total - supported

    warning_pipelines = sum(
        bool(result["warnings"])
        for result in pipeline_results
    )

    robustness = supported / total

    return {
        "supported": supported,
        "rejected": rejected,
        "total": total,
        "robustness": robustness,
        "warning_pipelines": warning_pipelines
    }
# ============================================================
# Output
# ============================================================
def save_results_csv(pipeline_results, filename):
    """
    Save pipeline-level results to a CSV file.

    Parameters
    ----------
    pipeline_results : list
        Results from all phylogenetic pipelines.

    filename : str
        Output CSV filename.
    """

    with open(filename, "w", newline="") as file:

        writer = csv.writer(file)

        # Header
        writer.writerow([
            "pipeline",
            "distance_model",
            "tree_method",
            "target_group",
            "conclusion",
            "warnings"
        ])

        # One row per pipeline
        for result in pipeline_results:

            target_group = ";".join(
                sorted(result["target_group"])
            )

            warnings = " | ".join(
                result["warnings"]
            )

            writer.writerow([
                result["pipeline"],
                result["distance_model"],
                result["tree_method"],
                target_group,
                result["conclusion"],
                warnings
            ])

# ============================================================
# Main program
# ============================================================

def main():

    args = parse_arguments()

    fasta_path = args.fasta

    # --------------------------------------------------
    # 1. Read input sequences
    # --------------------------------------------------

    sequences = read_fasta("Ver_01/test.fasta")

    # --------------------------------------------------
    # 2. Validate the alignment
    # --------------------------------------------------

    validate_alignment(sequences)
    os.makedirs("Ver_01/results", exist_ok=True)

    # --------------------------------------------------
    # 3. Define the biological claim being tested
    # --------------------------------------------------
     
    target_group = set(args.target)

    validate_target_group(
    sequences,
    target_group
)

    # --------------------------------------------------
    # 4. Run all four phylogenetic pipelines
    # --------------------------------------------------

    pipeline_results = []

    for pipeline in PIPELINES:

        result = analyse_pipeline(
            sequences,
            pipeline,
            target_group
        )

        pipeline_results.append(result)


# --------------------------------------------------
# Save Newick trees
# --------------------------------------------------


    for result in pipeline_results:

        filename = result["filename"] + ".nwk"

        filepath = os.path.join(
            "Ver_01",
            "results",
            filename
        )

        save_newick(
            result["tree"],
            filepath
    )
        print(f"Saved Newick tree: {filepath}")

# --------------------------------------------------
# Save pipeline results to CSV
# --------------------------------------------------

    csv_filepath = os.path.join(
        "Ver_01",
        "results",
        "pipeline_results.csv"
    )

    save_results_csv(
        pipeline_results,
        csv_filepath
    )

    print(f"Saved pipeline results: {csv_filepath}")
            

    

    # --------------------------------------------------
    # 5. Display individual pipeline results
    # --------------------------------------------------

    print("\nPipeline Results:")
    print("-" * 70)

    for result in pipeline_results:

        print(f"Pipeline: {result['pipeline']}")
        print(f"Distance model: {result['distance_model']}")
        print(f"Tree method: {result['tree_method']}")
        print(f"Target group: {result['target_group']}")
        print(f"Conclusion: {result['conclusion']}")

        if result["warnings"]:
            print("Warnings:")
            for warning in result["warnings"]:
                print(f"  - {warning}")
        else:
            print("Warnings: None")

        print("-" * 70)

    # --------------------------------------------------
    # 6. Calculate overall robustness
    # --------------------------------------------------

    robustness_result = calculate_robustness(
    pipeline_results)


    # --------------------------------------------------
    # 7. Display robustness
    # --------------------------------------------------

    print("\nRobustness Summary:")
    print("-" * 70)

    print(f"Target group: {target_group}")
    print(
        f"Pipelines tested: "
        f"{robustness_result['total']}"
    )

    print(
        f"Supported: "
        f"{robustness_result['supported']}"
    )

    print(
        f"Rejected: "
        f"{robustness_result['rejected']}"
    )

    print(
        f"Robustness: "
        f"{robustness_result['robustness']:.2%}"
    )

    print(
        f"Pipelines with warnings: "
        f"{robustness_result['warning_pipelines']}"
    )

    print("-" * 70)

if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        print(f"Error: {error}")