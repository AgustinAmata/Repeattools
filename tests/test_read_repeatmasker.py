import unittest
from pathlib import Path

from src.input_reader import read_repeatmasker

class RepMask(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_repeatmasker.txt"

    def test_read_repmask(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_repeatmasker(input_fhand)
            
        assert read_repeats[0] == {
            "sw": "1306", "per div": "15.6", "per del": "6.2",
            "per ins": "0.0", "query": "HSU08988", "q begin": "6563",
            "q end": "6781", "q left": "(22462)", "repeat": "C",
            "class/family": "MER7A", "r begin": "DNA/MER2_type",
            "r end": "(0)", "r left": "336", "id": "103"
            }
        
if __name__ == "__main__":
    unittest.main()