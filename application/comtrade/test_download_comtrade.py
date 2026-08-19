import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


MODULE_PATH = Path(__file__).with_name("download_comtrade.py")
SPEC = importlib.util.spec_from_file_location("download_comtrade", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
download = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = download
SPEC.loader.exec_module(download)


class SpecificationTests(unittest.TestCase):
    def test_requested_dimensions(self):
        self.assertEqual(len(download.COUNTRIES), 25)
        self.assertEqual(len({country.code for country in download.COUNTRIES}), 25)
        self.assertEqual(len(download.HS2_CODES), 96)
        self.assertNotIn("77", download.HS2_CODES)
        self.assertEqual(download.HS2_CODES[0], "01")
        self.assertEqual(download.HS2_CODES[-1], "97")

    def test_france_uses_comtrade_area_code(self):
        france = next(country for country in download.COUNTRIES if country.iso3 == "FRA")
        self.assertEqual(france.code, 251)

    def test_current_comtrade_customs_area_codes(self):
        codes = {country.iso3: country.code for country in download.COUNTRIES}
        self.assertEqual(codes["CHE"], 757)
        self.assertEqual(codes["IND"], 699)

    def test_full_plan(self):
        periods = download.month_range("2012-01", "2024-12")
        queries = download.build_queries("2012-01", "2024-12")
        self.assertEqual(len(periods), 156)
        self.assertEqual(len(queries), 325)
        self.assertTrue(all(len(query.periods) == 12 for query in queries))
        self.assertTrue(all(query.expected_cell_bound == 28_800 for query in queries))
        self.assertTrue(all(query.expected_cell_bound < download.MAX_RECORDS for query in queries))

    def test_sdk_settings_match_ui_specification(self):
        query = download.build_queries("2024-01", "2024-01")[0]
        settings = download.sdk_kwargs(query)
        self.assertEqual(settings["typeCode"], "C")
        self.assertEqual(settings["freqCode"], "M")
        self.assertEqual(settings["clCode"], "HS")
        self.assertEqual(settings["flowCode"], "X")
        self.assertEqual(settings["partner2Code"], "0")
        self.assertEqual(settings["customsCode"], "C00")
        self.assertEqual(settings["motCode"], "0")
        self.assertEqual(settings["breakdownMode"], "classic")
        self.assertIsNone(settings["aggregateBy"])
        self.assertEqual(len(settings["partnerCode"].split(",")), 25)
        self.assertEqual(len(settings["cmdCode"].split(",")), 96)

    def test_invalid_periods_fail_early(self):
        with self.assertRaises(ValueError):
            download.month_range("2024-13", "2024-12")
        with self.assertRaises(ValueError):
            download.month_range("2024-02", "2024-01")
        with self.assertRaises(ValueError):
            list(download.batched(("202401",), 13))

    def test_combine_uses_manifest_and_keeps_primary_value(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            relative = Path("reporter=842/reporter-842_USA_periods-202401-202401.csv")
            source = root / relative
            source.parent.mkdir(parents=True)
            source.write_text(
                "period,reporterCode,partnerCode,cmdCode,primaryValue\n"
                "202401,842,124,01,12345\n",
                encoding="utf-8",
            )
            source.with_suffix(".json").write_text(
                json.dumps({"status": "complete", "row_count": 1}), encoding="utf-8"
            )
            (root / "manifest.json").write_text(
                json.dumps({"checkpoint_csvs": [str(relative)]}), encoding="utf-8"
            )
            args = SimpleNamespace(
                output_dir=str(root),
                combined_file="combined.csv",
                overwrite=False,
                allow_partial=False,
            )
            self.assertEqual(download.run_combine(args), 0)
            combined = (root / "combined.csv").read_text(encoding="utf-8")
            self.assertIn("primaryValue", combined)
            self.assertIn("12345", combined)

    def test_unavailable_checkpoint_is_resumable(self):
        with tempfile.TemporaryDirectory() as temporary:
            query = next(
                query
                for query in download.build_queries("2012-01", "2012-12")
                if query.reporter.iso3 == "KOR"
            )
            root = Path(temporary)
            metadata_path = download.write_unavailable(root, query)
            csv_path, expected_metadata_path = download.completion_paths(root, query)
            self.assertEqual(metadata_path, expected_metadata_path)
            self.assertFalse(csv_path.exists())
            self.assertTrue(download.is_complete(csv_path, metadata_path, query))
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            self.assertEqual(metadata["status"], "unavailable")
            self.assertEqual(len(metadata["unavailable_periods"]), 12)

    def test_legacy_checkpoint_is_refreshed_for_partial_availability(self):
        with tempfile.TemporaryDirectory() as temporary:
            query = download.build_queries("2012-01", "2012-12")[0]
            root = Path(temporary)
            csv_path, metadata_path = download.completion_paths(root, query)
            csv_path.parent.mkdir(parents=True)
            csv_path.write_text("primaryValue\n1\n", encoding="utf-8")
            metadata_path.write_text(
                json.dumps(
                    {
                        "status": "complete",
                        "reporter_code": query.reporter.code,
                        "periods": list(query.periods),
                        "row_count": 1,
                        "parameters": download.sdk_kwargs(query),
                    }
                ),
                encoding="utf-8",
            )
            self.assertTrue(
                download.is_complete(
                    csv_path, metadata_path, query, available_periods=query.periods
                )
            )
            self.assertFalse(
                download.is_complete(
                    csv_path, metadata_path, query, available_periods=query.periods[3:]
                )
            )


if __name__ == "__main__":
    unittest.main()
