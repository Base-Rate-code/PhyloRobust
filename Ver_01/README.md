# PhyloRobust

**PhyloRobust V1.0** is a lightweight Python tool for assessing the robustness of a phylogenetic conclusion across alternative distance models and tree-building methods.

Instead of relying on a single phylogenetic pipeline, PhyloRobust asks:

> **How stable is a biological conclusion under reasonable methodological changes?**

---

## V1.0 Scope

PhyloRobust V1.0 evaluates a user-defined target clade using four analytical pipelines:

| Distance model | Tree method |
|---|---|
| p-distance | UPGMA |
| p-distance | Neighbor Joining |
| Jukes-Cantor | UPGMA |
| Jukes-Cantor | Neighbor Joining |

The target conclusion is considered **robust** when the target clade is supported across the tested pipelines.

---

## Features

- Aligned nucleotide FASTA input
- FASTA validation
- Target-group validation
- p-distance calculation
- Jukes-Cantor distance calculation
- Pairwise distance matrix generation
- UPGMA tree construction
- Neighbor Joining tree construction
- Target-clade detection
- Negative branch-length warnings
- Newick tree export
- Pipeline-level CSV results
- Overall robustness calculation
- Automated test suite

---

## Requirements

- Python 3
- Standard Python library only

No external Python packages are required.

---

## Usage

From the repository root:

```powershell
python Ver_01/phylo_robust_v1.py Ver_01/tests/test.fasta --target Sequence_A Sequence_B