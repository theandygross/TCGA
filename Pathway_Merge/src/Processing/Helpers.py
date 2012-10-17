'''
Created on Jul 12, 2012

@author: agross
'''
import pandas as pandas
from pandas import notnull
from numpy.linalg import LinAlgError
from numpy import diag

transferIndex = lambda source,target: pandas.Series(list(target), 
                                                    index=source.index)

from rpy2.robjects import r, FloatVector
def bhCorrection(s): 
    return pandas.Series(r('p.adjust')(FloatVector(s), method='BH'), 
                         index=s.index, name=s.name)
    
def match_series(a,b):
    a, b = a.align(b, join='inner', copy=False)
    valid = notnull(a) & notnull(b)
    a = a[valid]
    b = b[valid]
    return a,b

def split_a_by_b(a,b):
    a, b = match_series(a, b)
    groups = [a[b==num] for num in set(b)]
    return groups

def frame_svd(data_frame):
    '''
    Wrapper for taking in a pandas DataFrame, preforming SVD
    and outputting the U, S, and vH matricies in DataFrame form.
    '''
    from numpy.linalg import svd #@UnresolvedImport
    from pandas import DataFrame
    
    U,S,vH = svd(data_frame.as_matrix(), full_matrices=False)
    U = DataFrame(U, index=data_frame.index)
    vH = DataFrame(vH, columns=data_frame.columns).T
    return U,S,vH

def extract_pc(data_frame, pc_threshold=.2):
    try:
        U,S,vH = frame_svd(((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T)
    except LinAlgError:
        return None
    p = S**2/sum(S**2)
    return vH[0] if p[0] > pc_threshold else None

def drop_first_norm_pc(data_frame):
    '''
    Normalize the data_frame by rows and then reconstruct it without the first 
    principal component.  (Idea is to drop the biggest global pattern.)
    '''
    norm = ((data_frame.T - data_frame.mean(1)) / data_frame.std(1)).T
    U,S,vH = frame_svd(norm.dropna())
    S[0] = 0   #zero out first pc
    rest = U.dot(pandas.DataFrame(diag(S)).dot(vH.T))
    return rest