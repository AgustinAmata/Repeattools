import ast
import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from src.read_input import merge_inputs

class MergeIn(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path1 = test_path / "test_merge_inputs_rm.csv"
        self.test_path2 = test_path / "test_merge_inputs_te.csv"

    def test_mer_in(self):
        with open(self.test_path1) as rm_input_fhand, open(self.test_path2) as te_input_fhand:
            rm_input_df = pd.read_csv(rm_input_fhand)
            te_input_df = pd.read_csv(te_input_fhand, converters={"domains": ast.literal_eval})
            merged_df = merge_inputs(rm_input_df, te_input_df)

        convert_dict = {
        "class": "category", "superfamily": "category",
        "repeat": "category", "tes order":"category",
        "tes superfamily": "category", "clade": "category"
        }
        test_df= pd.DataFrame([{
            "per div": 30.6,"per del": 1.4,
            "per ins": 1.4, "seqid": "Peame105C00", "start": 10027900,
            "end": 10028180,
            "repeat": "rnd-5_family-987", "class": "LINE",
            "superfamily": "L1", "length": 280, 
            "tes order": "Unknown", "tes superfamily": "Unknown", "clade": "Unknown",
            "domains": [{"none": "none"}]
            }, {
            "per div": 30.6, "per del": 1.4,
            "per ins": 1.4, "seqid": "Peame105C00", "start": 10027969,
            "end": 10028180,
            "repeat": "rnd-5_family-987", "class": "LINE",
            "superfamily": "L1", "length": 211, "tes order": "LINE",
            "tes superfamily": "Unknown", "clade": "LINE",
            "domains": [{"RT": "LINE"}]
            }]).astype(convert_dict)
        
        assert_frame_equal(merged_df, test_df)

        
if __name__ == "__main__":
    unittest.main()