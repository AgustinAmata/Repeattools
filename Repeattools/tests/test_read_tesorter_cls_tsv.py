import unittest
from csv import DictReader
from pathlib import Path

from src.read_input import read_tesorter_cls_tsv

class TESort(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_tesorter_cls_tsv.tsv"

    def test_read_tes_tsv(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_tesorter_cls_tsv(input_fhand)
        assert read_repeats[0] == {
            "seqid": "Peame105C00", "start": "10027969", "end": "10028180",
            "repeat": "rnd-5_family-987", "class": "LINE",
            "superfamily": "L1", "tes order": "LINE",
            "tes superfamily": "unknown", "clade": "unknown",
            "complete": "unknown", "strand": "+", "domains": {"RT": "LINE"}
            }
        
if __name__ == "__main__":
    unittest.main()