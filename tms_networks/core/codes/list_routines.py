import numpy as np

def flattenLL(LL):
    ''' Flattens a list of lists: each element of each list
        becomes an element of a single 1-D list. '''
    return [x for L0 in LL for x in L0]

def nanpercentage(x):
    ''' Percentage of nan values from all values of
        the x array/N-dimensional matrix. '''
    return (np.sum(np.isnan(x))) / float(x.size)

