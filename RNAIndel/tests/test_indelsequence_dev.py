#!/usr/bin/env python3 

from unittest import TestCase

import sys
sys.path.append('..')
from  indelsequence_dev import Indel

class TestIndel(TestCase):
    
    def setUp(self):
        self.indel_1 = Indel('X', 200, 1, 'A')
        self.indel_2 = Indel('chrX', 200, 1, 'A')

    def test_chr(self):
        self.assertEqual(self.indel_1.chr, 'chrX')
        self.assertEqual(self.indel_2.chr, 'chrX')

class 


if __name__ == '__main__':
    from  unittest import main
    main()
