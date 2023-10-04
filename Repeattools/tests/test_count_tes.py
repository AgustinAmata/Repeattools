import unittest
from pathlib import Path

import pandas as pd
from pandas.testing import assert_series_equal

from src.create_matrix import count_tes

class CountTransposableElement(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_count_tes.csv"

    def test_count_tes(self):
        with open(self.test_path) as input_fhand:
            input_fhand = pd.read_csv(input_fhand)
            counted_tes = count_tes(input_fhand, "Persea_americana")

        test_series = pd.Series({"L1": 3}, name="Persea_americana", dtype="int32")
        test_series.index.name = "superfamily"

        assert_series_equal(counted_tes, test_series) 
        
if __name__ == "__main__":
    unittest.main()