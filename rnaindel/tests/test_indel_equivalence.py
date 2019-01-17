#!/usr/bin/env python3

from unittest import TestCase

try:
    from rnaindel.rnaindel_lib import SequenceWithIndel
except:
    from ..rnaindel_lib import SequenceWithIndel


class TestIndelEquivalentSolver(TestCase):
    
    def setUp(self):
        # insertion equivalent case 1
        self.idl1 = SequenceWithIndel('chr9', 37002671, 1, 'TTGGTG', 'TCGGCCTCCTACCCC', 'TCGGCG')
        self.idl2 = SequenceWithIndel('chr9', 37002676, 1, 'GTCGGC', 'CTCCTACCCCTCGGC', 'GCTGGG')
        # insertion equivalent case 2
        self.idl3 = SequenceWithIndel('chr22', 23653976, 1, 'GCGTGT', 'CCGG', 'CCGGTG')
        self.idl4 = SequenceWithIndel('chr22', 23653980, 1, 'GTCCGG', 'CCGG', 'TGTGGC')
        # insertion equivalent case 3
        self.idl5 = SequenceWithIndel('chr20', 34240741, 1, 'GGGGCA', 'GGGCCG', 'GGGCCGGGGCCGGGGCC')
        self.idl6 = SequenceWithIndel('chr20', 34240758, 1, 'GGGCCGGGGCCGGGGCC', 'GGGGCC', 'GGGGCC')

        
        # insertion equivalent case 3
        #self.idl5 = SequenceWithIndel('chr20', 34240741, 1, 'TCG', 'GGCTCAGGCTCA', 'GGCTCAGGCTCAGGCTC')
        #self.idl6 = SequenceWithIndel('chr20', 34240758, 1, 'GGGCTCAGGCTCAGGCTC', 'AGGCTCAGGCTC', 'TGCTTCGTACACTGG')

        
        
        # deletion equivalent case 1
        self.idl7 = SequenceWithIndel('chr4', 103827697, 0, 'ACCAAC', 'AAT', 'ATTTGA')
        self.idl8 = SequenceWithIndel('chr4', 103827698, 0, 'CCAACA', 'ATA', 'TTGAT')
        # deletion equivalent case 2
        self.idl9 = SequenceWithIndel('chrX', 140994892, 0, 'CAGAGC', 'CCT', 'CCTCAG')
        self.idl10 = SequenceWithIndel('chrX', 140994895, 0, 'AGCCCT', 'CCT', 'CAGGGG')
        # deletion equivalent case 3
        self.idl11 = SequenceWithIndel('chr1', 6638843, 0, 'GCTGGCAGCTAACAC', 'GCT', 'GCTGCTGCTGCTGCTGCTTG')
        self.idl12 = SequenceWithIndel('chr1', 6638860, 0, 'GCTGCTGCTGCTGCTGC', 'TGC', 'TTGGGACTGCTGGCCTGT')
        
    def test_equivalence(self):
        self.assertEqual(self.idl1 == self.idl2, True)
        self.assertEqual(self.idl3 == self.idl4, True)
        self.assertEqual(self.idl5 == self.idl6, True)
        self.assertEqual(self.idl7 == self.idl8, True)
        self.assertEqual(self.idl9 == self.idl10, True)
        self.assertEqual(self.idl11 ==self.idl12, True)

if __name__ == '__main__':
    from unittest import main
    main()
