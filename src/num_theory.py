from functools import reduce
from time_wrap import timed

def crt(res_mod):
    """finds number given residues modulo mutually coprime moduli"""
    e = 0
    
    mod_prod = reduce(lambda a, b: a * b[1], res_mod, 1)
    
    for res, mod in res_mod:        
        no_mod = mod_prod // mod

        u, v = ext_euclidean(mod, no_mod)
        
        e += res * v * no_mod
        
    return e % mod_prod

def ext_euclidean(a, b):
    
    x, y = 0, 1
    u, v = 1, 0
    
    while a != 0:
        q, r = b // a, b % a
        m, n = x - u * q, y - v * q
        
        b, a = a, r
        x, y = u, v  
        u, v = m, n
        
    return x, y

def bin_ext_euclidean(a, b):
    """extended euclidean algorithm for the second parameter a power of two"""
    
    u = 0
    v = 1
    
    b >>= 1
    s = b
    
    while s > 0:
        
        s >>= 1
        
        if (v & 1) ^ 1:
            u >>= 1
            v >>= 1
        else:
            u = (u >> 1) + b
            v = (v + a) >> 1
            
    return u, v

def gcd(a, b):
    
    u, v =  ext_euclidean(a, b)
    
    return u * a + v * b

def is_coprime(a, b):
    
    return gcd(a, b) == 1

def is_mutually_coprime(a, b, c):
    
    return is_coprime(a, b) and is_coprime(b, c) and is_coprime(c, a)

def is_coprime_all(a, m_coprimes):
    
    for num in m_coprimes:
        if not is_coprime(a, num):
            return False
    
    return True

def odd_num_evens(a, b, c):
    
    ne = a % 2 + b % 2 + c % 2
    
    return ne % 2 == 0

#@timed
def exp_bs(a, p, m):
    """raises a to the pth power modulo m using exponentiation by squares"""
    
    a %= m
    r = 1
    
    while p:
        if p & 1:
            r = (r * a) % m
        a = (a * a) % m
        p >>= 1

    return r

def dlog_mod(g, h, m, mt, pe):
    
    p = pe[0]
    e = pe[1]
    
    f = p**e
    no_tot = mt // f
    
    aeg = exp_bs(g, no_tot, m)
    aeh = exp_bs(h, no_tot, m)
    oeg = exp_bs(g, mt // p, m)
    
    r = 0
    d = 1
    c = 0
    
    while f != 1:
        f //= p
        
        ntg = exp_bs(aeg, f, m)
        nth = exp_bs(aeh, f, m)
        
        q = exp_bs(ntg, r, m)
        
        for i in range(p):
            if q == nth:
                break
            else:
                q = q * oeg % m
        
        r += i * d                
        d *= p
        c += 1
        
        #print("\tx = {!s} (mod {!s}^{!s})".format(r, p, c))
        
    return r

def dlog(g, h, m, mf):
    """finds discrete log of rh base rg modulo m using Pohlig-Hellman algorithm"""
    
    x_crt = []
    
    for pe in mf[2]:
        r = dlog_mod(g, h, m, mf[0], pe)
        x_crt.append([r, pe[0]**pe[1]])
        
        #print("x = {!s} (mod {!s}^{!s})".format(r, pe[0], pe[1]))
    
    return crt(x_crt)

#techincally don't need all of mf, but included for compatability with dlog    
def dlog_brute(g, h, m, mf):      
    
    e = mf[0]
    
    l = 1
    
    for i in range(e):            
        if l == h:
            break
        else:
            l = l * g % m
    else:
        return None
    
    return i    

def is_gen_n_k(g, m, mf):
    """returns whether g is a generator of the multiplicative group modulo n^k for all positive k"""
    
    if not is_coprime(g, m):
        return False
    
    tot = mf[0]
    fact = mf[1]
    tot_fact = mf[2]
    
    #determine whether g a generator modulo n^1
    for p, e in tot_fact:
        if exp_bs(g, tot // p, m) == 1:
            return False        
    
    #determine whether g a generator modulo p^k for all p | n and k > 1
    #this implies g a generator modulo n^k if g a generator modulo n^1
    for p, e in fact:
        sm = p**2
        
        if exp_bs(g, p - 1, sm) == 1:
            return False
    
    return True

@timed
def p_count_gen_n_k(m):
    """prints a count of the number of generators modulo m**k for all positive k"""
    
    c = 0
    
    for g in range(m):
        if is_gen_n_k(g, m):
            c += 1
    
    print("There are {!s} generators modulo {!s}**k for all positive k".format(c, m))

def p_test_gen_n_k(g, m, mf, km):
    """brute force tests g modulo m**k for k=1, ..., km and compares the number of unique values to the totient function"""
    
    tot = mf[0]
    
    print(is_gen_n_k(g, m))
    
    r = 1
    
    for k in range(1, km + 1):        
        cp_set = set()
        
        r *= m
        l = 1
    
        for i in range(tot):
            l %= r
            cp_set.add(l)
            l *= g
        
        print("Expected {!s} coprime members modulo {!s}**{!s}".format(tot, m, k))
        print("Received {!s} coprime members modulo {!s}**{!s}".format(len(cp_set), m, k))
        
        tot *= m
