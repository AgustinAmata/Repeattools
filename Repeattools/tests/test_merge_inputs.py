import ast
import unittest
from pathlib import Path

from src.read_input import merge_inputs

class MergeIn(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path1 = test_path / "test_merge_inputs_rm.txt"
        self.test_path2 = test_path / "test_merge_inputs_te.txt"

    def test_mer_in(self):
        with open(self.test_path1) as rm_input_fhand, open(self.test_path2) as te_input_fhand:
            rm_input_fhand = ast.literal_eval(rm_input_fhand.read())
            te_input_fhand = ast.literal_eval(te_input_fhand.read())
            merged_inputs = merge_inputs(rm_input_fhand, te_input_fhand)
        assert merged_inputs[0] == {
            "sw": "452", "per div": "30.6","per del": "1.4",
            "per ins": "1.4", "seqid": "Peame105C00", "start": "10027900",
            "end": "10028180", "q left": "(45826121)", "match": "+",
            "repeat": "rnd-5_family-987", "class": "LINE",
            "superfamily": "L1", "r start": "659", "r end": "870",
            "r left": "(467)", "id": "7765", "domains": None,
            "tes order": None, "tes superfamily": None, "complete": None,
            "strand": None, "clade": None
            }
        assert merged_inputs[1] == {
            "sw": "452", "per div": "30.6", "per del": "1.4",
            "per ins": "1.4", "seqid": "Peame105C00", "start": "10027969",
            "end": "10028180", "q left": "(45826121)", "match": "+",
            "repeat": "rnd-5_family-987", "class": "LINE",
            "superfamily": "L1", "r start": "659", "r end": "870",
            "r left": "(467)", "id": "7765", "tes order": "LINE",
            "tes superfamily": "unknown", "clade": "unknown",
            "complete": "unknown", "strand": "+", "domains": {"RT": "LINE"}
            }

        
if __name__ == "__main__":
    unittest.main()