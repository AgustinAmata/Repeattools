import ast
import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from src.create_matrix import create_te_count_matrix

class CreateMatrix(unittest.TestCase):

    def test_create_matrix(self):
        input_fhand1 = pd.Series({"L1": 3}, name="Peame104")
        input_fhand1.index.name = "superfamily"
        input_fhand2 = pd.Series({"SINE": 3}, name="Peame103")
        input_fhand2.index.name = "superfamily"

        te_count_matrix = create_te_count_matrix([input_fhand1, input_fhand2])

        test_df = pd.DataFrame({"Peame104": [3.0, 0.0], "Peame103": [0.0, 3.0]}, index=["L1", "SINE"])

        assert_frame_equal(te_count_matrix, test_df) 
        
if __name__ == "__main__":
    unittest.main()