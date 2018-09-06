#!/usr/bin/env python3

from unittest import TestCase

import sys
sys.path.append('../rna_indel_lib')
from most_common import most_common

class MostCommon(TestCase):

   def test_most_common(self):
       self.assertRaises(ValueError, most_common, [])
       self.assertRaises(ValueError, most_common, None)
       self.assertEqual(most_common(['a', 'c', 'd', 'd', 'a', 'd', 'c']), 'd')
       self.assertEqual(most_common(['b', 'a', 'c', 'a', 'c']) in ['a', 'c'], True)  
                 
if __name__ == '__main__':
    from unittest import main
    main()
