import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from src.read_input import read_repeatmasker_out

class RepMask(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_repeatmasker_out.out"

    def test_read_repmask(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_repeatmasker_out(input_fhand)
        convert_dict = {
        "per div": "float16", "per del": "float16",
        "per ins": "float16","start": "int32", "end": "int32",
        "repeat": "category", "length": "int32", "class": "category",
        "superfamily": "category"
        }
        test_df = pd.DataFrame([{
            "per div": "14.7", "per del": "2.1", "per ins": "2.4",
            "seqid": "Peame105C00", "start": "23", "end": "4959",
            "repeat": "rnd-4_family-1935", "length": 4936, "class": "Unknown",
            "superfamily": "Unknown"
            },{
            "per div": "31.9", "per del": "5.0", "per ins": "5.0",
            "seqid": "Peame105C00", "start": "4959", "end": "5519",
            "repeat": "ltr-1_family-724", "length": 560, "class": "LTR",
            "superfamily": "Copia"
            },{
            "per div": "32.2", "per del": "4.6", "per ins": "4.6",
            "seqid": "Peame105C00", "start": "5630", "end": "5847",
            "repeat": "ltr-1_family-645",  "length": 217, "class": "LTR",
            "superfamily": "Copia"
            }]).astype(convert_dict)
        
        assert_frame_equal(test_df, read_repeats)
if __name__ == "__main__":
    unittest.main()