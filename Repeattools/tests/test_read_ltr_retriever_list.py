import unittest
from csv import DictReader
from pathlib import Path

from src.input_reader import read_ltr_retriever_list

class LTRRetList(unittest.TestCase):
    
    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_ltr_retriever_list.txt"

    def test_ltr_ret_ls(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_ltr_retriever_list(input_fhand)

        assert read_repeats[0] == {
            "#LTR_loc": "Chr10:100..200", "Category": "pass",
            "Motif":"motif:TGCC", "TSD": "TSD:GGTGG",
            "5'_TSD":"94..99", "3'_TSD": "201..206",
            "Internal": "IN:105..195", "Similarity": "0.9875",
            "Strand": "-", "Family": "Gypsy", "Superfamily": "LTR",
            "Insertion_Time": "10000"
            }

if __name__ == "__main__":
    unittest.main()