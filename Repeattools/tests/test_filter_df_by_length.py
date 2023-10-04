import ast
import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from src.create_matrix import filter_df_by_length

class FilterLength(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_filter_df_by_length.csv"

    def test_filter_length(self):
        with open(self.test_path) as input_fhand:
            input_fhand = pd.read_csv(input_fhand, converters={"domains": ast.literal_eval})
            filtered_df = filter_df_by_length(input_fhand, 220)

        test_df = pd.DataFrame([
            {"per div": 30.6, "per del": 1.4,"per ins": 1.4,
             "seqid": "Peame105C00", "start": 10027900,"end": 10028180,
             "repeat": "rnd-5_family-987",
             "class": "LINE","superfamily": "L1",
             "domains": [{"none": "none"}], "tes order": "Unknown",
             "tes superfamily": "Unknown",
             "clade": "Unknown", "length": 280},
             ])

        assert_frame_equal(filtered_df.reset_index(drop=True), test_df.reset_index(drop=True)) 


        
if __name__ == "__main__":
    unittest.main()