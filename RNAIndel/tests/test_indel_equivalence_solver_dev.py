#!/usr/bin/env python3

from unittest import TestCase

import sys
sys.path.append('..')
import rna_indel_lib as rna


class TestIndelEquivalentSolver(TestCase):
    
    def setUp(self):
        # insertion equivalent case 1
        self.idl1 = rna.SequenceWithIndel('chr9', 37002671, 1, 'TTGGTG', 'TCGGCCTCCTACCCC', 'TCGGCG')
        self.idl2 = rna.SequenceWithIndel('chr9', 37002676, 1, 'GTCGGC', 'CTCCTACCCCTCGGC', 'GCTGGG')
        # insertion equivalent case 2
        self.idl3 = rna.SequenceWithIndel('chr22', 23653976, 1, 'GCGTGT', 'CCGG', 'CCGGTG')
        self.idl4 = rna.SequenceWithIndel('chr22', 23653980, 1, 'GTCCGG', 'CCGG', 'TGTGGC')
        # insertion equivalent case 3
        self.idl5 = rna.SequenceWithIndel('chr20', 34240741, 1, 'GGGGCA', 'GGGCCG', 'GGGCCGGGGCCGGGGCC')
        self.idl6 = rna.SequenceWithIndel('chr20', 34240758, 1, 'GGGCCGGGGCCGGGGCC', 'GGGGCC', 'GGGGCC')

        # deletion equivalent case 1
        self.idl7 = rna.SequenceWithIndel('chr4', 103827697, 0, 'ACCAAC', 'AAT', 'ATTTGA')
        self.idl8 = rna.SequenceWithIndel('chr4', 103827698, 0, 'CCAACA', 'ATA', 'TTGAT')
        # deletion equivalent case 2
        self.idl9 = rna.SequenceWithIndel('chrX', 140994892, 0, 'CAGAGC', 'CCT', 'CCTCAG')
        self.idl10 = rna.SequenceWithIndel('chrX', 140994895, 0, 'AGCCCT', 'CCT', 'CAGGGG')
        # deletion equivalent case 3
        self.idl11 = rna.SequenceWithIndel('chr1', 6638843, 0, 'GCTGGCAGCTAACAC', 'GCT', 'GCTGCTGCTGCTGCTGCTTG')
        self.idl12 = rna.SequenceWithIndel('chr1', 6638860, 0, 'GCTGCTGCTGCTGCTGC', 'TGC', 'TTGGGACTGCTGGCCTGT')
        
    def test_are_equivalent(self):
        self.assertEqual(rna.are_equivalent(self.idl1, self.idl2), True)
        self.assertEqual(rna.are_equivalent(self.idl3, self.idl4), True)
        self.assertEqual(rna.are_equivalent(self.idl5, self.idl6), True)
        self.assertEqual(rna.are_equivalent(self.idl7, self.idl8), True)
        self.assertEqual(rna.are_equivalent(self.idl9, self.idl10), True)
        self.assertEqual(rna.are_equivalent(self.idl11, self.idl12), True)

if __name__ == '__main__':
    from unittest import main
    main()
