import unittest
from csv import DictReader
from pathlib import Path

from src.read_input import read_ltr_retriever_list

class LTRRetList(unittest.TestCase):
    
    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_ltr_retriever_list.txt"

    def test_ltr_ret_ls(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_ltr_retriever_list(input_fhand)

        assert read_repeats[0] == {
            "seqid": "Chr10", "start": "100", "end": "200", "category": "pass",
            "motif":"motif:TGCC", "tsd": "TSD:GGTGG",
            "5'_tsd":"94..99", "3'_tsd": "201..206",
            "internal": "IN:105..195", "similarity": "0.9875",
            "strand": "-", "family": "Gypsy", "superfamily": "LTR",
            "insertion_time": "10000"
            }

if __name__ == "__main__":
    unittest.main()