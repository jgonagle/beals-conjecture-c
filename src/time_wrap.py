from time import time
from functools import wraps
import ipdb

def timed(func):
    
    @wraps(func)
    def time_this(*args, **kargs):
        
        s = time()
        
        r = func(*args, **kargs)
        
        t = time() - s
        
        print("{} took {!s} seconds".format(func.__name__, t))
        
        return r
    
    return time_this