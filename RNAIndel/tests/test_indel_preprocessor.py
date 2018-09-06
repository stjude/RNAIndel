#!/usr/bin/env python3

import pandas as pd
from unittest import TestCase
from unittest.mock import patch
from pandas.util.testing import assert_frame_equal

import sys
sys.path.append('../rna_indel_lib')
import indel_preprocessor as pre


class TestIndelPreprocessor(TestCase):
    
    @patch('indel_preprocessor.exists_bambino_output')
    def test_exists_bambino_output(self, mock_exists):
        
        # mocking the case where bambino output exists
        mock_exists.return_value = True
        pre.exists_bambino_output('bambino_output.txt')
        self.assertEqual(pre.exists_bambino_output('bambino.txt'), True)
        
        # mocking the case where bambino output does not exist
        mock_exists.return_value = False
        pre.exists_bambino_output('no_such_file.txt')
        self.assertEqual(pre.exists_bambino_output('no_such_file.txt'), False)
     

    def test_is_canonical_chromosome(self):
        self.assertEqual(pre.is_canonical_chromosome(pd.Series({'Chr': 'chr6_cox_hap1'})), False)
        self.assertEqual(pre.is_canonical_chromosome(pd.Series({'Chr': 'chrY'})), True)
        self.assertEqual(pre.is_canonical_chromosome(pd.Series({'Chr': 'chr3'})), True)

    def setUp(self):
        df = pd.read_csv('./testdata/expected_1.txt', sep='\t')
        self.fixture = df
     
    def test_main(self):
        assert_frame_equal(pre.indel_preprocessor('./testdata/test_1.txt'), self.fixture)
               
if __name__ == '__main__':
    from unittest import main
    main()
