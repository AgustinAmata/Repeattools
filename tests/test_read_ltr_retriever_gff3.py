import unittest
from csv import DictReader
from pathlib import Path

from src.input_reader import read_ltr_retriever_gff3

class LTRRetgff3(unittest.TestCase):

    def setUp(self):
        test_path = Path(__file__).parent.absolute()
        test_path = test_path / "data"
        self.test_path = test_path / "test_read_ltr_retriever_gff3.txt"

    def test_read_ltr_ret_gff3(self):
        with open(self.test_path) as input_fhand:
            read_repeats = read_ltr_retriever_gff3(input_fhand)

        assert read_repeats[0] == {
            "seqid": "Tzi8_chr1", "source": "LTR_retriever",
            "repeat_class/superfamily": "repeat_region",
            "start": "54481", "end": "63973", "sw_score": ".",
            "strand": "-", "phase": ".",
            "attributes": "ID=repeat_region1;Sequence_ontology=SO:0000657;Name=Tzi8_chr1:54486..63968#LTR/Copia;motif=TGCA;tsd=GAGTG;ltr_identity=1.0000;Method=structural"
            }

if __name__ == "__main__":
    unittest.main()