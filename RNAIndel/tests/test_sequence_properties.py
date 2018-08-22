#!/usr/bin/env python3

from  unittest import TestCase

import sys
sys.path.append('..')
import sequence_properties as sp


class SequenceProperties(TestCase):
    """ test for sequence_properties.py"""
    
    def test_editdistance(self):
        self.assertEqual(sp.editdistance('', ''), 0) 
        self.assertEqual(sp.editdistance('AAT', ''), 3)
        self.assertEqual(sp.editdistance('GAAGCA', 'AAGCAA'), 2)

    def test_linguistic_complexity(self):
        self.assertRaises(ValueError, sp.linguistic_complexity, 'NAGG')
        self.assertEqual(sp.linguistic_complexity(''), 0)
        # 'GCTGCT'
        #  i = 1, max_vocab = min(4, 6-1+1) = 4 -> 3/4
        #  i = 2, max_vocab = min(4^2, 6-2+1) = 5 -> 3/5
        #  i = 3, max_vocab = min(4^3, 6-3+1) = 4 -> 3/4
        #  i = 4, max_vocab = min(4^4, 6-4+1) = 3 -> 3/3
        #  i = 5, max_vocab = min(4^5, 6-5+1) = 2 -> 2/2
        #  i = 6, max_vocab = 1 -> always 1 (no need to calculate) 
        self.assertAlmostEqual(sp.linguistic_complexity('GCTGCT'), (3*3*3*3*2)/(4*5*4*3*2))

    def test_repeat(self):
        self.assertRaises(ValueError, sp.repeat, None, 'AAT', 'T', 'TTTTTG')
        self.assertRaises(ValueError, sp.repeat, 2, 'AAT', 'T', 'TTTTTG')
        self.assertRaises(ValueError, sp.repeat, 1, '', 'A', 'TACG')
        self.assertRaises(ValueError, sp.repeat, 1, 'GT', 'A', None)
        self.assertEqual(sp.repeat(0, 'GTAGAG', 'AGAGAGAG', 'AGAGTC'), 8)
        self.assertEqual(sp.repeat(1, 'GTAGAG', 'AGAGAGAG', 'AGAGTC'), 4)
        self.assertEqual(sp.repeat(0, 'GTAGAG', 'AGAGAGA', 'AGAGTC'), 0)
    
    def test_dna_strength(self):
        self.assertRaises(ValueError, sp.dna_strength, None)
        self.assertRaises(ValueError, sp.dna_strength, 'G')
        # 'AGCGTGA'
        # AG = 8, GC = 13, CG = 10, GT = 10, TG = 7, GA = 8
        self.assertEqual(sp.dna_strength('AGCGTGA'), (8+13+10+10+7+8)/len('AGCGTGA'))
    
    def test_gc(self):
        self.assertRaises(ValueError, sp.gc, None)
        self.assertRaises(ValueError, sp.gc, '')
        self.assertEqual(sp.gc('TTTAAA'), 0)
        self.assertEqual(sp.gc('ATGATC'), 2/6)

    def test_dissimilarity(self):
        self.assertRaises(ValueError, sp.dissimilarity, 'ATGAC', '', 'GTAT')
        self.assertRaises(ValueError, sp.dissimilarity, None, 'AACTG', 'GTAT')
        self.assertEqual(sp.dissimilarity('T', 'ATGCA', 'ATGCA'), 0)
        self.assertEqual(sp.dissimilarity('GGAGGA', 'TGT', 'CCTCGA'), 2/len('TGT'))

    def test_reverse_complement(self):
        self.assertRaises(ValueError, sp.reverse_complement, None)
        self.assertRaises(ValueError, sp.reverse_complement, '')
        self.assertEqual(sp.reverse_complement('NNNN'), 'NNNN')
        self.assertEqual(sp.reverse_complement('ATCG'), 'CGAT')

    def test_exists_stop_codon(self):
        self.assertRaises(ValueError, sp.exists_stop_codon, None, 'ATGACT')
        self.assertRaises(ValueError, sp.exists_stop_codon, 1, 'ATGACT')
        self.assertEqual(sp.exists_stop_codon('-', 'CAGTTACAT'), True) 
        self.assertEqual(sp.exists_stop_codon('+', 'CAGTTAAAT'), False)

if __name__ == '__main__':
    from unittest import main
    main()   
