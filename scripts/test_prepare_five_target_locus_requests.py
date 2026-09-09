#!/usr/bin/env python3

import unittest

from scripts import prepare_five_target_locus_requests as prep


class FiveTargetLocusRequestTests(unittest.TestCase):
    def test_panel_has_expected_factor_and_matrix_counts(self) -> None:
        tracks = prep.score_tracks()
        self.assertEqual(len(tracks), 30)
        self.assertEqual(len({track["track_id"] for track in tracks}), 30)
        self.assertEqual(len({track["factors"][0]["factor_id"] for track in tracks}), 28)
        tfap2c = [track for track in tracks if track["factors"][0]["factor_id"] == "TFAP2C"]
        self.assertEqual([track["source_ids"][0] for track in tfap2c],
                         ["MA0524.3", "MA0814.3", "MA0815.1"])
        self.assertEqual([track["source_ids"][0] for track in tracks if
                          track["factors"][0]["factor_id"] == "KLF15"], ["MA1513.2"])

    def test_target_panel_is_deduplicated(self) -> None:
        self.assertEqual(set(prep.TARGET_TRANSCRIPTS), {"CD44", "TGFB1", "SERPINE1", "PATZ1", "TP73"})
        self.assertTrue(all(len(ids) == len(set(ids)) == 3 for ids in prep.TARGET_TRANSCRIPTS.values()))


if __name__ == "__main__":
    unittest.main()
