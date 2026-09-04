import unittest
import sys
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1])
)

import phylo_robust_v1 as pr

class TestDistanceModels(unittest.TestCase):

    def test_p_distance_identical(self):
        self.assertEqual(
            pr.p_distance("AAAA", "AAAA"),
            0.0
        )

    def test_p_distance_known_difference(self):
        self.assertEqual(
            pr.p_distance("AAAA", "AAAT"),
            0.25
        )
    def test_jukes_cantor_identical(self):
        self.assertAlmostEqual(
            pr.jukes_cantor_distance(
                "AAAA",
                "AAAA"
            ),
            0.0
        )

    def test_jukes_cantor_known_difference(self):
        expected = -0.75 * __import__("math").log(
            1 - (4 / 3) * 0.25
        )

        self.assertAlmostEqual(
            pr.jukes_cantor_distance(
                "AAAA",
                "AAAT"
            ),
            expected
        )

class TestValidation(unittest.TestCase):

    def test_valid_alignment(self):
        sequences = {
            "A": "AAAA",
            "B": "AAAT",
            "C": "AATT"
        }

        self.assertTrue(
            pr.validate_alignment(sequences)
        )

    def test_unequal_lengths(self):
        sequences = {
            "A": "AAAA",
            "B": "AAA"
        }

        with self.assertRaises(ValueError):
            pr.validate_alignment(sequences)

    def test_valid_target(self):
        sequences = {
        "A": "AAAA",
        "B": "AAAT",
        "C": "AATT"
    }

        self.assertIsNone(
            pr.validate_target_group(
                sequences,
                {"A", "B"}
            )
        )

    def test_missing_target(self):
        sequences = {
            "A": "AAAA",
            "B": "AAAT",
            "C": "AATT"
        }

        with self.assertRaises(ValueError):
            pr.validate_target_group(
                sequences,
                {"A", "X"}
            )

class TestRobustness(unittest.TestCase):

    def test_all_supported(self):

        results = [
            {"target_supported": True, "warnings": []},
            {"target_supported": True, "warnings": []},
            {"target_supported": True, "warnings": []},
            {"target_supported": True, "warnings": []}
        ]

        result = pr.calculate_robustness(results)

        self.assertEqual(result["supported"], 4)
        self.assertEqual(result["rejected"], 0)
        self.assertEqual(result["total"], 4)
        self.assertEqual(result["robustness"], 1.0)
        self.assertEqual(result["warning_pipelines"], 0)

    def test_mixed_results(self):

        results = [
            {"target_supported": True, "warnings": []},
            {"target_supported": True, "warnings": ["warning"]},
            {"target_supported": False, "warnings": []},
            {"target_supported": True, "warnings": []}
        ]

        result = pr.calculate_robustness(results)

        self.assertEqual(result["supported"], 3)
        self.assertEqual(result["rejected"], 1)
        self.assertEqual(result["total"], 4)
        self.assertEqual(result["robustness"], 0.75)
        self.assertEqual(result["warning_pipelines"], 1)

class TestPipelineAnalysis(unittest.TestCase):

    def setUp(self):
        self.sequences = {
            "Sequence_A": "AAAA",
            "Sequence_B": "AAAT",
            "Sequence_C": "AATT"
        }

        self.target_group = {
            "Sequence_A",
            "Sequence_B"
        }

    def test_all_four_pipelines(self):

        for pipeline in pr.PIPELINES:

            result = pr.analyse_pipeline(
                self.sequences,
                pipeline,
                self.target_group
            )

            self.assertIn(
                "pipeline",
                result
            )

            self.assertIn(
                "tree",
                result
            )

            self.assertIn(
                "target_supported",
                result
            )

            self.assertIn(
                "warnings",
                result
            )

            self.assertTrue(
                result["target_supported"]
            )

    def test_jukes_cantor_nj_warning(self):

        pipeline = next(
            p for p in pr.PIPELINES
            if p["name"] == "Jukes-Cantor + NJ"
        )

        result = pr.analyse_pipeline(
            self.sequences,
            pipeline,
            self.target_group
        )

        self.assertTrue(
            result["target_supported"]
        )

        self.assertGreater(
            len(result["warnings"]),
            0
        )

        self.assertIn(
            "Negative branch length detected",
            result["warnings"][0]
        )   

class TestDistanceMatrix(unittest.TestCase):

    def setUp(self):
        self.sequences = {
            "Sequence_A": "AAAA",
            "Sequence_B": "AAAT",
            "Sequence_C": "AATT"
        }

    def test_p_distance_matrix(self):

        matrix = pr.calculate_distance_matrix(
            self.sequences,
            pr.p_distance
        )

        self.assertAlmostEqual(
            matrix["Sequence_A"]["Sequence_B"],
            0.25
        )

        self.assertAlmostEqual(
            matrix["Sequence_A"]["Sequence_C"],
            0.50
        )

        self.assertAlmostEqual(
            matrix["Sequence_B"]["Sequence_C"],
            0.25
        )

    def test_distance_matrix_symmetry(self):

        matrix = pr.calculate_distance_matrix(
            self.sequences,
            pr.p_distance
        )

        self.assertEqual(
            matrix["Sequence_A"]["Sequence_B"],
            matrix["Sequence_B"]["Sequence_A"]
        )

        self.assertEqual(
            matrix["Sequence_A"]["Sequence_C"],
            matrix["Sequence_C"]["Sequence_A"]
        )

class TestNewick(unittest.TestCase):

    def test_target_clade_in_newick(self):

        sequences = {
            "Sequence_A": "AAAA",
            "Sequence_B": "AAAT",
            "Sequence_C": "AATT"
        }

        pipeline = next(
            p for p in pr.PIPELINES
            if p["name"] == "p-distance + UPGMA"
        )

        result = pr.analyse_pipeline(
            sequences,
            pipeline,
            {"Sequence_A", "Sequence_B"}
        )

        newick = pr.export_newick(
            result["tree"]
        )

        self.assertTrue(
            newick.endswith(";")
        )

        self.assertIn(
        "(Sequence_A:0.125,Sequence_B:0.125)",
        newick)
    

        self.assertIn(
        "Sequence_C:0.1875",
        newick)

        self.assertEqual(
         newick,
        "(Sequence_C:0.1875,(Sequence_A:0.125,Sequence_B:0.125):0.0625);"
        )
class TestEdgeCases(unittest.TestCase):

    def test_empty_target(self):

        sequences = {
            "A": "AAAA",
            "B": "AAAT"
        }

        with self.assertRaises(ValueError):
            pr.validate_target_group(
                sequences,
                set()
            )

    def test_missing_target_taxon(self):

        sequences = {
            "A": "AAAA",
            "B": "AAAT"
        }

        with self.assertRaises(ValueError):
            pr.validate_target_group(
                sequences,
                {"A", "X"}
            )

    def test_target_single_taxon(self):

        sequences = {
            "A": "AAAA",
            "B": "AAAT",
            "C": "AATT"
        }

        # A single taxon is not a meaningful clade target
        with self.assertRaises(ValueError):
            pr.validate_target_group(
                sequences,
                {"A"}
            )    

    def test_empty_alignment(self):

        sequences = {}

        with self.assertRaises(ValueError):
            pr.validate_alignment(sequences)   

        def test_invalid_nucleotide(self):

            sequences = {
                "A": "AAAA",
                "B": "AAAX"
            }

            with self.assertRaises(ValueError):
                pr.validate_alignment(sequences)

if __name__ == "__main__":
    unittest.main()