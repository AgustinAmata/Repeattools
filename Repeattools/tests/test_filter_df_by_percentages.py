import ast
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from src.create_matrix import filter_df_by_percentages

class FilterPerc(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_filter_df_by_percentages.csv"

    def test_filter_perc(self):
        with open(self.test_path) as input_fhand:
            input_fhand = pd.read_csv(input_fhand, converters={"domains": ast.literal_eval})
            filtered_df = filter_df_by_percentages(input_fhand, "ins=4.0")

        test_df = pd.DataFrame([
            {"sw": "452", "per div": "30.6","per del": "1.4","per ins": "1.4",
             "seqid": "Peame105C00", "start": "10027900","end": "10028180",
             "q left": "(45826121)", "match": "+","repeat": "rnd-5_family-987",
             "class": "LINE","superfamily": "L1", "r start": "659",
             "r end": "870","r left": "(467)", "id": "7765",
             "domains": [{"none": "none"}], "tes order": "none",
             "tes superfamily": "none", "complete": "none","strand": "none",
             "clade": "none"},
              {"sw": "452", "per div": "30.6", "per del": "1.4",
              "per ins": "3.0", "seqid": "Peame105C00", "start": "10027969",
              "end": "10028180", "q left": "(45826121)", "match": "+",
              "repeat": "rnd-5_family-987", "class": "LINE",
              "superfamily": "L1", "r start": "659", "r end": "870",
              "r left": "(467)", "id": "7765", "tes order": "LINE",
              "tes superfamily": "unknown", "clade": "Ale",
              "complete": "unknown", "strand": "+", 
              "domains": [{"RT": "Ale"}, {"RH": "Ale"}]}
              ])
        
        convert_dict = {
        "sw": "int64", "per div": "float64", "per del": "float64",
        "per ins": "float64","start": "int64", "end": "int64", "r start": "int64",
        "r end": "int64", "id": "int64"
        }
        test_df = test_df.astype(convert_dict)

        assert_frame_equal(filtered_df.reset_index(drop=True), test_df.reset_index(drop=True)) 
        
if __name__ == "__main__":
    unittest.main()