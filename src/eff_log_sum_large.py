from math import log, exp
from time import time
from time_wrap import timed
import decimal as d
import ipdb

precision = 100
o = d.Decimal(1)
d.getcontext().prec = precision

def is_between_large(a, x, b, y, c, z):
    
    e = x * precise_log(a)
    f = y * precise_log(b)
    g = z * precise_log(c)
    
    return ibq_large(e, f, g)

def ibq_large(e, f, g):
    
    lb, ub = eff_lsq_large(e, f)
    
    if lb <= g <= ub:
        return True
    
    return False

#@timed
def eff_log_sum_large(a, x, b, y):
    
    e = x * precise_log(a)
    f = y * precise_log(b)
    
    return eff_lsq_large(e, f)

def eff_lsq_large(e, f):
    
    g, h = d.Decimal(e).max(f), d.Decimal(e).min(f)
    
    r = precise_exp(h - g)
    
    lb = g + r * (o - r * (o / 2 - r * (o / 3 - r * (o / 4 - r * (o / 5 - r / 6)))))
    ub = g + r * (o - r * (o / 2 - r * (o / 3 - r * (o / 4 - r / 5))))
    
    return [lb, ub]

#@timed
def normal_log_sum_large(a, x, b, y):
    
    r = a**x + b**y
    e = precise_log(r)
    
    return e

def precise_log(x):
    
    return d.Decimal(x).ln()

def precise_exp(x):
    
    return d.Decimal(x).exp()

def precise_logminus(max, min):
    
    if max <= min:
        return 0
        
    return d.Decimal(max) + (1 - d.Decimal(min - max).exp()).ln()

def precise_logplus(max, min):
            
    return d.Decimal(max) + (1 + d.Decimal(min - max).exp()).ln()

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
                            
                                if ibq_large(a, x, b, y, c, z):
                                    print("{!s}^{!s} + {!s}^{!s} ~= {!s}^{!s}".format(a, x, b, y, c, z))
                    