import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from src.read_input import read_tesorter_cls_tsv

class TESort(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_tesorter_cls_tsv.tsv"

    def test_read_tes_tsv(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_tesorter_cls_tsv(input_fhand)

        convert_dict = {
        "clade": "category", "tes superfamily": "category",
        "tes order": "category", "start": "int32", "end": "int32",
        "repeat": "category", "length": "int32", "class": "category",
        "superfamily": "category"
        }
        test_df = pd.DataFrame([{
            "tes order": "LINE", "tes superfamily": "Unknown",
            "clade": "Unknown", "domains": [{"RT": "LINE"}],
            "seqid": "Peame105C00", "start": "10027969", "end": "10028180",
            "repeat": "rnd-5_family-987", "class": "LINE",
            "superfamily": "L1", "length": 211
            },{
            "tes order": "LTR", "tes superfamily": "Copia",
            "clade": "Ale", "domains": [{"GAG": "Ale"}],
            "seqid": "Peame105C00", "start": "10037928", "end": "10038407",
            "repeat": "ltr-1_family-331", "class": "LTR",
            "superfamily": "Copia", "length": 479
            },{
            "tes order": "LTR", "tes superfamily": "Copia",
            "clade": "Ale", "domains": [{"RT": "Ale"},{"RH": "Ale"}],
            "seqid": "Peame105C00", "start": "10038812", "end": "10040267",
            "repeat": "ltr-1_family-437", "class": "LTR",
            "superfamily": "Unknown", "length": 1455
            }]).astype(convert_dict)
        
        assert_frame_equal(test_df, read_repeats)

if __name__ == "__main__":
    unittest.main()