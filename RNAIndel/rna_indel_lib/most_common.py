#!/usr/bin/env python3


def most_common(lst):
    """Returns the most common element in a list

    Args:
        lst (list): list with any datatype
    Returns:
        an element of lst (any): the most common element
                                 if commonests tie, randomly selected 
    Raises:
       ValueError: if lst is None or empty list
    """
    if not lst:
        raise ValueError("input must be non-empty list")

    return max(set(lst), key=lst.count)
