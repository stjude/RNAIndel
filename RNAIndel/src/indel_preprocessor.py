#!/usr/bin/env python3

import os
import sys
import logging
import pandas as pd

logger = logging.getLogger(__name__)

def main(filename):
    """Validate, extract and format indel calls from Bambino output
     
    Args:
        filename (str): Bambino output filename (contains SNVs + indels)
    Returns:
        df (pandas.DataFrame): Only contains indels
                               4 columns formatted as: 

                               chr  pos  ref  alt
                        
                               chr1 123  -    ATC         
                               chr2 456  GG   -
                                      ....
                               chrY 987  CCT  -
    """
    if not exists_bambino_output(filename):
        sys.exit(1)
    
    df = pd.read_csv(filename, sep='\t', dtype={'dbSNP': str})
     
    if not is_non_trivial_result(df):
        sys.exit(1)
    
    df = extract_necessary_info(df)
    df = extract_indel_calls(df)
    
    if len(df) == 0:
        logging.warning('No indels detected in variant calling. Analysis done.')
        sys.exit(0)
        
    df = rename_header(df)
    df = format_indel_report(df)
    df = df.reset_index(drop=True)
    
    return df


def exists_bambino_output(filename):
    """Assert if Bambino output file exists
       
    Args:
        filename (str): Bambino output filename
    Returns:
        it_exists (bool): True if exists
    """
    it_exists = False
           
    if not os.path.exists(filename):
        logging.critical('Bambino output does NOT exist.')
    else:
        it_exists = True

    return it_exists

    
def is_non_trivial_result(df):
    """Assert if Bambino outputs contains data other than the header line

    Args:
        df (pandas.DataFrame): Bambino output as pd.DataFrame
    Returns:
        is_non_trivial (bool): True if it has more than header line
    """
    is_non_trivial = False

    if len(df) == 0:
        logging.critical('Bambino output only contains the header line.')
    else:
        is_non_trivial = True

    return is_non_trivial


def extract_necessary_info(df):
    """Reduce data size by selecting columns and chromosomes

    Args:
        df (pandas.DataFrame): 45 columns
    Retuns:
        df (pandas.DataFrame): 5 columns with canonical chromosomes
    """
    pd.options.mode.chained_assignment = None

    df = df[['Chr', 'Pos', 'Type', 'Chr_Allele', 'Alternative_Allele']]
    
    df['is_canonical'] = df.apply(is_canonical_chromosome, axis=1)
    df = df[df['is_canonical'] == 1]
    df.drop('is_canonical', axis=1, inplace=True)

    return df 


def is_canonical_chromosome(row):
    """Check if chr is 1-22, X or Y (M not included)

    Args:
        row (pandas.Series): chromosome name (str) at index 'Chr'
    Returns:
        is_canonical (bool): True if chr is 1-22, X or Y
    """
    is_canonical = False
    
    chr_name = row['Chr'].replace('chr', '')
    
    if chr_name == 'X' or chr_name == 'Y':
        is_canonical = True
    else:
        try:
            chr_num = int(chr_name)
            if 1 <= chr_num <= 22:
                is_canonical = True
        except:
            pass
    
    return is_canonical


def extract_indel_calls(df):
    """Extract indels calls and sort by position

    Args: 
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df['original_order'] = df.index

    df_d = df[df['Type'] == 'deletion']
    df_i = df[df['Type'] == 'insertion']
    df = pd.concat([df_d, df_i])
    
    df.sort_values('original_order', inplace=True)
    df.drop('original_order', axis=1, inplace=True)
     
    return df 


def rename_header(df):
    """Rename as follows
       Chr -> chr
       Pos -> pos
       Chr_Allele -> ref
       Alternative_Allele -> alt

    Args:
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df.drop('Type', axis=1, inplace=True) 
    df = df.rename(columns={'Chr': 'chr',
                            'Pos': 'pos',
                            'Chr_Allele': 'ref',
                            'Alternative_Allele': 'alt'})
    
    return df

                            
def format_indel_report(df):
    """Format as follows
    For insertion
           ref          alt
           -            inserted_seq
    For deletion
           ref          alt
           deleted_seq  -
    Args:
        df (pandas.DataFrame)
    Returns:
        df (pandas.DataFrame)
    """
    df['ref'] = df.apply(lambda x: '-' if x['ref'] != x['ref'] else x['ref'], axis=1)
    df['alt'] = df.apply(lambda x: '-' if x['alt'] != x['alt'] else x['alt'], axis=1) 
    
    return df
