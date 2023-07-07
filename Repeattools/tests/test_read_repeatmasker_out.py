import unittest
from pathlib import Path

from src.read_input import read_repeatmasker_out

class RepMask(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_repeatmasker_out.out"

    def test_read_repmask(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_repeatmasker_out(input_fhand)
        assert read_repeats[0] == {
            "sw": "25785", "per div": "14.7", "per del": "2.1",
            "per ins": "2.4", "seqid": "Peame105C00", "start": "23",
            "end": "4959", "q left": "(55849342)", "match": "+",
            "repeat": "rnd-4_family-1935", "class": "Unknown",
            "superfamily": "Unknown", "r start": "1", "r end": "4924",
            "r left": "(0)", "id": "1"
            }
        
if __name__ == "__main__":
    unittest.main()