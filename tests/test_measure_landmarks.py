"""Synthetic arithmetic and file-safety checks; no clinical validation."""

import copy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "measure_landmarks.py"
SPEC = importlib.util.spec_from_file_location("measure_landmarks", SCRIPT)
calculator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(calculator)


def record(**changes):
    value = {
        "id": "SYN-TEST",
        "kind": "length",
        "pixel_spacing_mm": [2, 0.5],
        "ipp_mm": [10, 20, 30],
        "iop": [1, 0, 0, 0, 1, 0],
        "points": [[0, 0], [6, 2]],
    }
    value.update(changes)
    return value


def document(*records):
    return {"coordinate_convention": "zero_based_pixel_centers",
            "measurements": list(records) or [record()]}


class GeometryTests(unittest.TestCase):
    def assertVectorAlmostEqual(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for result, target in zip(actual, expected):
            self.assertAlmostEqual(result, target)

    def assertFiniteTree(self, value):
        if isinstance(value, dict):
            for child in value.values():
                self.assertFiniteTree(child)
        elif isinstance(value, (tuple, list)):
            for child in value:
                self.assertFiniteTree(child)
        elif isinstance(value, (float, int)):
            self.assertTrue(math.isfinite(value))

    def test_published_synthetic_fixture(self):
        source = json.loads((ROOT / "references" / "synthetic-landmarks.json").read_text())
        results = calculator.calculate_document(source)["measurements"]
        self.assertEqual([r["computed_value_mm"] for r in results], [3, 5, 4])
        self.assertFiniteTree(results)

    def test_anisotropic_spacing_and_pixel_center_origin(self):
        result = calculator.calculate(record())
        self.assertEqual(result["computed_value_mm"], 5)
        self.assertEqual(result["landmarks_patient_coordinates_mm"],
                         [(10, 20, 30), (13, 24, 30)])

    def test_oblique_orientation_preserves_physical_length(self):
        diagonal = math.sqrt(0.5)
        result = calculator.calculate(record(
            iop=[diagonal, diagonal, 0, -diagonal, diagonal, 0]))
        self.assertAlmostEqual(result["computed_value_mm"], 5)
        self.assertVectorAlmostEqual(result["landmarks_patient_coordinates_mm"][1],
                                     [10 - diagonal, 20 + 7 * diagonal, 30])

    def test_fractional_centers_and_negative_patient_origin(self):
        result = calculator.calculate(record(
            points=[[0, 0], [0.5, 0.5]], rows=2, columns=2,
            ipp_mm=[-10, -20, -30]))
        self.assertEqual(result["landmarks_patient_coordinates_mm"][1],
                         (-9.75, -19, -30))
        self.assertAlmostEqual(result["computed_value_mm"], math.hypot(0.25, 1))

    def test_axis_projection_preserves_signed_direction(self):
        for points, expected in (([[0, 0], [6, 2]], -4), ([[6, 2], [0, 0]], 4)):
            with self.subTest(points=points):
                result = calculator.calculate(record(
                    kind="patient_axis_projection", points=points,
                    iop=[1, 0, 0, 0, 0, -1], axis=[0, 0, 1]))
                self.assertEqual(result["signed_axis_difference_mm"], expected)
                self.assertEqual(result["computed_value_mm"], 4)

    def test_oblique_axis_projection_is_not_segment_length(self):
        diagonal = math.sqrt(0.5)
        result = calculator.calculate(record(
            kind="patient_axis_projection", axis=[diagonal, diagonal, 0]))
        self.assertAlmostEqual(result["computed_value_mm"], 7 * diagonal)

    def test_perpendicular_is_distance_to_infinite_line(self):
        for row, parameter, outside in ((0, -0.2, True), (1, 0, False),
                                         (4, 0.6, False), (6, 1, False),
                                         (8, 1.4, True)):
            with self.subTest(row=row):
                result = calculator.calculate(record(
                    kind="perpendicular", line=[[1, 1], [1, 6]], point=[7, row]))
                self.assertEqual(result["computed_value_mm"], 3)
                self.assertAlmostEqual(result["projection_parameter"], parameter)
                self.assertEqual(result["foot_outside_reference_segment"], outside)
                self.assertVectorAlmostEqual(result["perpendicular_foot_patient_mm"],
                                             [10.5, 20 + 2 * row, 30])

    def test_large_finite_perpendicular_regression(self):
        source = record(kind="perpendicular", pixel_spacing_mm=[1, 1],
                        ipp_mm=[0, 0, 0], line=[[0, 0], [1e200, 0]],
                        point=[1e200, 0])
        result = calculator.calculate(source)
        self.assertEqual(result["computed_value_mm"], 0)
        self.assertEqual(result["projection_parameter"], 1)
        self.assertEqual(result["perpendicular_foot_patient_mm"], (1e200, 0, 0))
        self.assertFiniteTree(result)

    def test_nonfinite_derived_values_are_rejected_by_python_api(self):
        cases = {
            "scaled_offset": record(pixel_spacing_mm=[1, 1e308],
                                    points=[[0, 0], [2, 0]]),
            "patient_coordinate": record(pixel_spacing_mm=[1, 1],
                                         ipp_mm=[1e308, 0, 0],
                                         points=[[0, 0], [1e308, 0]]),
            "length": record(pixel_spacing_mm=[1, 1], ipp_mm=[0, 0, 0],
                             points=[[0, 0], [1.4e308, 1.4e308]]),
            "axis_projection": record(kind="patient_axis_projection",
                                      pixel_spacing_mm=[1, 1], ipp_mm=[0, 0, 0],
                                      points=[[0, 0], [1.4e308, 1.4e308]],
                                      axis=[math.sqrt(0.5), math.sqrt(0.5), 0]),
            "projection_parameter": record(kind="perpendicular",
                                           pixel_spacing_mm=[1, 1], ipp_mm=[0, 0, 0],
                                           line=[[0, 0], [1e-11, 0]], point=[1e308, 1]),
            "foot_coordinate": record(kind="perpendicular", pixel_spacing_mm=[1, 1],
                                      ipp_mm=[1e308, -1e308, 0],
                                      line=[[0, 0], [5e307, 5e307]], point=[0, 1.7e308]),
        }
        for label, source in cases.items():
            with self.subTest(stage=label):
                with self.assertRaisesRegex(ValueError, "non-finite|overflow"):
                    calculator.calculate(source)

    def test_numeric_inputs_reject_nonfinite_boolean_string_and_huge_integer(self):
        for invalid in (float("nan"), float("inf"), -float("inf"), True, "1", 10 ** 1000):
            for key, value in (("pixel_spacing_mm", [invalid, 1]),
                               ("ipp_mm", [invalid, 0, 0]),
                               ("iop", [invalid, 0, 0, 0, 1, 0]),
                               ("points", [[0, 0], [invalid, 1]])):
                with self.subTest(field=key, invalid_type=type(invalid).__name__):
                    with self.assertRaises(ValueError):
                        calculator.calculate(record(**{key: value}))

    def test_invalid_spacing_orientation_and_axis(self):
        invalid_changes = [
            {"pixel_spacing_mm": [0, 1]}, {"pixel_spacing_mm": [-1, 1]},
            {"iop": [2, 0, 0, 0, 1, 0]}, {"iop": [1, 0, 0, 1, 0, 0]},
            {"iop": [1, 0, 0]},
            {"kind": "patient_axis_projection", "axis": [0, 0, 0]},
            {"kind": "patient_axis_projection", "axis": [0, 0, 2]},
            {"kind": "patient_axis_projection", "axis": [0, 0, float("nan")]},
        ]
        for changes in invalid_changes:
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    calculator.calculate(record(**changes))

    def test_degenerate_reference_line(self):
        for end in ([0, 0], [1e-13, 0]):
            with self.subTest(end=end):
                with self.assertRaisesRegex(ValueError, "coincident"):
                    calculator.calculate(record(kind="perpendicular", ipp_mm=[0, 0, 0],
                                                line=[[0, 0], end], point=[3, 4]))

    def test_source_center_bounds(self):
        self.assertEqual(calculator.calculate(record(rows=3, columns=7))["computed_value_mm"], 5)
        for point in ([-0.01, 0], [0, -0.01], [6.01, 2], [6, 2.01]):
            with self.subTest(point=point):
                with self.assertRaises(ValueError):
                    calculator.calculate(record(rows=3, columns=7, points=[[0, 0], point]))

    def test_dimensions_must_be_a_pair_of_positive_integers(self):
        for changes in ({"rows": 3}, {"columns": 7}, {"rows": True, "columns": 7},
                        {"rows": 3.0, "columns": 7}, {"rows": 3, "columns": 0},
                        {"rows": -1, "columns": 7}):
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(ValueError, "positive integers"):
                    calculator.calculate(record(**changes))

    def test_optional_dimensions_do_not_claim_an_upper_bound(self):
        result = calculator.calculate(record(points=[[0, 0], [100, 0]]))
        self.assertEqual(result["computed_value_mm"], 50)

    def test_decimals_contract_applies_to_both_python_entry_points(self):
        for decimals in (-1, 5, True, 1.5, "2", None):
            with self.subTest(decimals=decimals):
                with self.assertRaisesRegex(ValueError, "decimals"):
                    calculator.calculate(record(), decimals)
                with self.assertRaisesRegex(ValueError, "decimals"):
                    calculator.calculate_document(document(), decimals)
        source = record(points=[[0, 0], [1, 1]])
        for decimals in range(5):
            result = calculator.calculate(source, decimals)
            self.assertEqual(result["display_mm"], round(math.hypot(0.5, 2), decimals))

    def test_preserves_synthetic_provenance_without_mutating_input(self):
        source = record(note="Synthetic provenance metadata")
        original = copy.deepcopy(source)
        result = calculator.calculate(source)
        self.assertEqual(source, original)
        self.assertEqual(result["input_record"], original)


class DocumentTests(unittest.TestCase):
    def test_rejects_missing_convention_empty_records_and_invalid_record(self):
        for value in ([], {}, {"measurements": [record()]},
                      {"coordinate_convention": "zero_based_pixel_centers", "measurements": []},
                      document({}), document(record(id=" ")), document(record(kind="volume"))):
            with self.subTest(value_type=type(value).__name__):
                with self.assertRaises(ValueError):
                    calculator.calculate_document(value)

    def test_duplicate_ids_are_scoped_to_study(self):
        with self.assertRaisesRegex(ValueError, "measurement 2: duplicate"):
            calculator.calculate_document(document(record(), record()))
        result = calculator.calculate_document(document(record(study="SYN-A"), record(study="SYN-B")))
        self.assertEqual(len(result["measurements"]), 2)
        with self.assertRaisesRegex(ValueError, "study"):
            calculator.calculate_document(document(record(study=1)))

    def test_error_identifies_failed_measurement_without_partial_result(self):
        with self.assertRaisesRegex(ValueError, "measurement 2:"):
            calculator.calculate_document(document(record(id="SYN-A"), record(id="SYN-B", points=[])))


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.workspace = tempfile.TemporaryDirectory(prefix="synthetic-landmarks-tests-")
        self.addCleanup(self.workspace.cleanup)
        self.directory = Path(self.workspace.name)
        self.input = self.directory / "synthetic-input.json"
        self.output = self.directory / "synthetic-output.json"
        self.input.write_text(json.dumps(document()), encoding="utf-8")

    def run_cli(self, *arguments):
        return subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), str(self.input), *arguments],
                              cwd=self.directory, env={"PATH": os.defpath},
                              capture_output=True, text=True, timeout=10)

    def assertNoOutputOrStagingFiles(self):
        self.assertFalse(self.output.exists())
        self.assertEqual(list(self.directory.glob(".landmarks-*.tmp")), [])

    def test_cli_stdout_is_complete_strict_json_with_input_hash(self):
        completed = self.run_cli()
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(completed.stdout)
        self.assertEqual(result["input_sha256"], hashlib.sha256(self.input.read_bytes()).hexdigest())
        self.assertEqual(result["measurements"][0]["computed_value_mm"], 5)
        self.assertNoOutputOrStagingFiles()

    def test_cli_creates_complete_private_output_and_does_not_overwrite_it(self):
        completed = self.run_cli("--output", str(self.output))
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(completed.stdout, "")
        previous = self.output.read_bytes()
        self.assertEqual(json.loads(previous)["measurements"][0]["computed_value_mm"], 5)
        if os.name == "posix":
            self.assertEqual(self.output.stat().st_mode & 0o777, 0o600)
        repeated = self.run_cli("--output", str(self.output))
        self.assertEqual(repeated.returncode, 2)
        self.assertEqual(self.output.read_bytes(), previous)
        self.assertEqual(list(self.directory.glob(".landmarks-*.tmp")), [])

    def test_cli_never_replaces_its_input_file(self):
        previous = self.input.read_bytes()
        completed = self.run_cli("--output", str(self.input))
        self.assertEqual(completed.returncode, 2)
        self.assertEqual(self.input.read_bytes(), previous)

    def test_cli_errors_leave_no_destination(self):
        cases = [
            "{invalid json",
            json.dumps(document(record(id="SYN-A"), record(id="SYN-B", points=[]))),
            json.dumps(document(record(note=float("nan")))),
            json.dumps(document(record(note="\ud800"))),
            json.dumps(document(record(pixel_spacing_mm=[1, 1e308], points=[[0, 0], [2, 0]]))),
        ]
        for index, content in enumerate(cases):
            with self.subTest(case=index):
                self.input.write_text(content, encoding="utf-8")
                completed = self.run_cli("--output", str(self.output))
                self.assertEqual(completed.returncode, 2)
                self.assertEqual(completed.stdout, "")
                self.assertNoOutputOrStagingFiles()

    def test_cli_rejects_invalid_decimals_before_creating_output(self):
        completed = self.run_cli("--output", str(self.output), "--decimals", "5")
        self.assertEqual(completed.returncode, 2)
        self.assertNoOutputOrStagingFiles()

    def test_existing_symlink_is_not_followed_or_replaced(self):
        self.output.symlink_to(self.input)
        previous = self.input.read_bytes()
        with self.assertRaises(FileExistsError):
            calculator.write_new_output(self.output, "new data")
        self.assertTrue(self.output.is_symlink())
        self.assertEqual(self.input.read_bytes(), previous)
        self.assertEqual(list(self.directory.glob(".landmarks-*.tmp")), [])

    def test_partial_staging_write_failure_leaves_no_destination(self):
        original = calculator.tempfile.NamedTemporaryFile

        def fail_after_some_bytes(**kwargs):
            stream = original(**kwargs)
            write = stream.write

            def partial_write(data):
                write(data[:5])
                raise OSError("synthetic disk-write failure")

            stream.write = partial_write
            return stream

        with mock.patch.object(calculator.tempfile, "NamedTemporaryFile", side_effect=fail_after_some_bytes):
            with self.assertRaisesRegex(OSError, "synthetic disk-write"):
                calculator.write_new_output(self.output, "longer synthetic output")
        self.assertNoOutputOrStagingFiles()

    def test_flush_failure_leaves_no_destination(self):
        with mock.patch.object(calculator.os, "fsync", side_effect=OSError("synthetic flush failure")):
            with self.assertRaisesRegex(OSError, "synthetic flush"):
                calculator.write_new_output(self.output, "complete synthetic payload")
        self.assertNoOutputOrStagingFiles()

    def test_competing_creation_is_preserved(self):
        def competing_creator(source, destination):
            destination.write_text("synthetic competing result", encoding="utf-8")
            raise FileExistsError("synthetic creation race")

        with mock.patch.object(calculator.os, "link", side_effect=competing_creator):
            with self.assertRaises(FileExistsError):
                calculator.write_new_output(self.output, "our synthetic result")
        self.assertEqual(self.output.read_text(), "synthetic competing result")
        self.assertEqual(list(self.directory.glob(".landmarks-*.tmp")), [])

    def test_unsupported_atomic_publish_leaves_no_destination(self):
        with mock.patch.object(calculator.os, "link", side_effect=OSError("synthetic unsupported hard link")):
            with self.assertRaisesRegex(OSError, "unsupported hard link"):
                calculator.write_new_output(self.output, "synthetic result")
        self.assertNoOutputOrStagingFiles()


if __name__ == "__main__":
    unittest.main()
