"""UI result regression: nested route evidence must not pass through CSV parsing."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from webui import server


class LargeEvidenceResults(unittest.TestCase):
    def test_preflight_reads_json_even_with_oversized_csv_field(self):
        for stage in ('risk', 'availability'):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                folder = root / stage
                folder.mkdir()
                rows = [{'id': 'generated-test', 'decision': 'review',
                         'route_evidence': {'detail': 'x' * 150000}}]
                (folder / 'records.json').write_text(json.dumps(rows), encoding='utf-8')
                with (folder / 'molecules.csv').open('w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(['id', 'route_evidence'])
                    writer.writerow(['generated-test', json.dumps(rows[0]['route_evidence'])])
                with patch.object(server.workflow, 'run_dir', return_value=root):
                    status, _, body = server.api_project_rows({}, 'test-run', stage)
                self.assertEqual(status, 200)
                self.assertEqual(json.loads(body)['rows'], rows)

    def test_preflight_json_does_not_depend_on_csv_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'availability').mkdir()
            rows = [{'id': 'candidate-a'}, {'id': 'candidate-b'}]
            (root / 'availability/records.json').write_text(json.dumps(rows))
            with patch.object(server.workflow, 'run_dir', return_value=root):
                _, _, body = server.api_project_rows({'search': ['candidate-b']}, 'test-run', 'availability')
            self.assertEqual(json.loads(body)['rows'], [rows[1]])

    def test_other_stage_csv_and_pagination_remain_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'reinvent4').mkdir()
            (root / 'reinvent4/molecules.csv').write_text('id,score\na,1\nb,2\n')
            with patch.object(server.workflow, 'run_dir', return_value=root):
                _, _, body = server.api_project_rows({'offset': ['1']}, 'test-run', 'reinvent4')
            self.assertEqual(json.loads(body), {'rows': [{'id': 'b', 'score': '2'}], 'total': 2, 'offset': 1})


if __name__ == '__main__':
    unittest.main()
