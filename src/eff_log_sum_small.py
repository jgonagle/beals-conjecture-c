from math import log, exp
from time import time
from time_wrap import timed
import ipdb

def is_between_small(a, x, b, y, c, z):
    
    e = x * log(a)
    f = y * log(b)
    g = z * log(c)
    
    return ibq_small(e, f, g)

def ibq_small(e, f, g):
    
    lb, ub = eff_lsq_small(e, f)
    
    if lb <= g <= ub:
        return True
    
    return False

#@timed
def eff_log_sum_small(a, x, b, y):
    
    e = x * log(a)
    f = y * log(b)
    
    return eff_lsq_small(e, f)

def eff_lsq_small(e, f):
    
    g, h = max(e, f), min(e, f)
    
    r = exp(h - g)
    
    lb = g + r * (1 - r * (1 / 2 - r * (1 / 3 - r * (1 / 4 - r * (1 / 5 - r / 6)))))
    ub = g + r * (1 - r * (1 / 2 - r * (1 / 3 - r * (1 / 4 - r / 5))))
    
    return [lb, ub]

#@timed
def normal_log_sum_small(a, x, b, y):
    
    r = a**x + b**y
    e = log(r)
    
    return e

#@timed
def normal_log_sum_approx(a, x, b, y):
    
    e = x * log(a)
    f = y * log(b)
    
    return max(e, f)

def test_all():
    
    min_power = 3
    log_two = log(2)
    
    for a in range(3, 100):
        print(a)
        j = log(a)
        
        for b in range(a, 100):
            print('\t' + str(b))
            k = log(b)
            
            for c in range(3, 100):
                print('\t\t' + str(c))
                l = log(c)
                e = j * min_power - 1
                
                if a != b and b != c and c != a:
                    for x in range(min_power, 100):
                        e += j
                        f = k * min_power - 1
                        
                        for y in range(min_power, 100):
                            f += k
                            g = l * min_power - 1
                            
                            for z in range(min_power, 100):
                                g += l
                            
                                if ibq_small(a, x, b, y, c, z):
                                    print("{!s}^{!s} + {!s}^{!s} ~= {!s}^{!s}".format(a, x, b, y, c, z))
                    