import math
import ipdb
from functools import reduce
from time_wrap import timed
from num_theory import crt, bin_ext_euclidean

"""all operations assumed on unsigned integers"""

class MontgomeryReduction:
    
    def __init__(self, m=None, tot_fact=None):
        
        if m is None:
            self.m = None
            self.t = None
            self.tf = None
                
            self.r = None
            self.mp = None
            self.rp = None
            self.l = None
            self.n = None
            
            self.ro = None
        else:
            self.set_modulus(m, tot_fact)
    
    def set_modulus(self, m, tot_fact=None):
        
        self.m = m
        
        if tot_fact is None:
            self.t = None
            self.tf = None
        else:
            self.t = tot_fact[0]
            self.tf = tot_fact[2]
            
        self.set_next_two_power()
        self.set_bin_inverse()
        
        self.ro = self.convert(1)
        
    #@timed
    def convert(self, x):
        """converts modulus form to montgomery reduced form
        
        x must not be montgomery reduced"""
        
        x <<= self.n
        x %= self.m
    
        return x
    
    #@timed
    def mult(self, a, b):
        """performs montgomery multiplication
        
        a and b must be montgomery reduced
        returns a result that is montgomery reduced"""
        
        t = a * b
        k = self.mult_mod_r(t, self.mp)
        
        u = (t + k * self.m) >> self.n
        u = u - (self.m & -(u >= self.m))
        
        return u
    
    #@timed
    def exp(self, a, p):
        """raises a to the pth power modulo m
        
        a must be montgomery reduced"""
        
        r = self.convert(1)
        
        while p:
            if p & 1:
                r = self.mult(r, a)
            a = self.mult(a, a)
            p >>= 1
    
        return r
    
    def dlog_mod(self, g, h, pe):
        
        p = pe[0]
        e = pe[1]
        
        f = p**e
        no_tot = self.t // f
        
        aeg = self.exp(g, no_tot)
        aeh = self.exp(h, no_tot)
        oeg = self.exp(g, self.t // p)
        
        r = 0
        d = 1
        c = 0
        
        while f != 1:
            f //= p
            
            ntg = self.exp(aeg, f)
            nth = self.exp(aeh, f)
            
            q = self.exp(ntg, r)
            
            for i in range(p):
                if q == nth:
                    break
                else:
                    q = self.mult(q, oeg)
            
            r += i * d                
            d *= p
            c += 1
            
            #print("\tx = {!s} (mod {!s}^{!s})".format(r, p, c))
            
        return r
    
    def dlog(self, g, h):
        """finds discrete log of rh base rg modulo m using Pohlig-Hellman algorithm"""
        
        if self.t is None:
            return None
        
        x_crt = []
        
        for pe in self.tf:
            r = self.dlog_mod(g, h, pe)
            x_crt.append([r, pe[0]**pe[1]])
            
            #print("x = {!s} (mod {!s}^{!s})".format(r, pe[0], pe[1]))
        
        return crt(x_crt)
    
    def dlog_brute(self, g, h):      
        
        e = (m - 1) if self.t is None else self.t
        
        l = self.ro
        
        for i in range(e):            
            if l == h:
                break
            else:
                l = self.mult(l, g)
        else:
            return None
        
        return i
    
    #@timed
    def inv(self, x):
        """returns x back to modulus form from montgomery reduced form
        
        x must be montgomery reduced"""
        
        return (x * self.rp) % self.m 
    
    def mult_mod_r(self, a, b):
        """multiply two numbers modulo r, a power of two"""
        
        a &= self.l
        b &= self.l
        
        s = (a * b) & self.l
        
        return s
    
    def set_next_two_power(self):
        
        l = self.m
        s = 1
        
        while -~l & l:
            l |= l >> s
            s *= 2
        
        r = l + 1
        
        n = 0
        t = l
        
        while t:
            t >>= 1
            n += 1
        
        self.r = r
        self.l = l
        self.n = n

    def set_bin_inverse(self):
        
        self.mp, self.rp = bin_ext_euclidean(self.m, self.r)