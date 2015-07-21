import math
import ipdb
from random import randrange
from time_wrap import timed
from populate_num_lists import get_all_lists
from num_theory import is_coprime, is_gen_n_k, crt
from montgomery_reduction import MontgomeryReduction

int_ceiling = 1000
max_pow = 10**3

num_ps_tests = 5

prob_primes = []
safe_primes = []
all_fact = []
prime_sq_fact = []

def try_all_beals(lim=1000):
    
    populate_all(lim)
    
    #b is in the primes since composite moduli don't have primitve roots.  Techincally, b could be in the prime powers or 2p^k, but the first is redundant for this application while the second isn't satisfied since b must be a power greater than two
    for b in prob_primes:
        for a in range(3, lim):
            if is_coprime(a, b) and is_gen_n_k(a, b, all_fact[b]):            
                for c in range(3, lim):
                    if odd_num_evens(a, b, c) and is_coprime(a, c) and is_coprime(b, c):
                        for z in range(3, lim):
                            r = beals_solver(a, b, c, z)
                            
                            print(a, b, c, z)
                            
                            if r[2][0] and r[0] > 2 and r[1] > 2:
                                print("{!s}^{!s} {} {!s}^{!s} = {!s}^{!s}".format(a, r[0], r[2][1], b, r[1], c, z))
                                wait = input("Possible solution found.  Press enter to continue.")

@timed
def test_beals_solver(n_tests=100):
    
    populate_all()
    
    num_pp = len(prob_primes)
    num_af = len(all_fact)
    
    print("Beginning {!s} tests of Beals Solver...".format(n_tests))
    
    while (n_tests > 0):
        
        a = randrange(num_af)
        b = prob_primes[randrange(num_pp)]
        
        if is_gen_n_k(a, b, all_fact[b]):
            x = randrange(num_af)
            
            #need to ensure that y is large enough that phi(b^y) > x
            min_y = math.ceil(math.log(num_af / (b - 1), b)) + 1
            y = randrange(min_y, num_af)
            
            ax = a**x
            by = b**y
            
            if randrange(0, 2):
                c = ax + by
                s = '+'
            else:
                c = ax - by
                
                if c < 0:
                    continue
                else:
                    s = '-'
            
            #print("\Solving {!s}^{!s} {!s} {!s}^{!s}...".format(a, x, s, b, y))
            r = beals_solver(a, b, c, 1)
            
            assert r[0] == x and r[1] == y and r[2][0] and r[2][1] == s
            
            n_tests -= 1
            
    print("All tests passed!")            

#@timed
def beals_solver(a, b, c, z):
    """finds x and y s.t. a^x +/- b^y = c^z
    
    a, b, c mutually coprime and a is a generator modulo b
    b should be prime, but we'll accept (very) strong probable primes to make finding generators modulo b easier to find, b will be a safe prime
    c must be positive"""
    
    x = 0
    y = 0
    
    #take ceiling for rounding errors (possible underflow)
    lrab = math.ceil(math.log(a, b))
    max_yc = z * math.ceil(math.log(c, b))
    
    b_fact = all_fact[b]
    #for modulus b^k for k > 1, we want the discrete log to be over a "totient" of b, since the possible values of x satisfying g^x=h (mod b^k) will range over 0 to b.  Giving the factorization of b as the factorization of the totient allows us to use the speedup afforded by the Pohlig-Hellman algorithm.  The second element is None because we don't need the factorization of the inverse totient of b
    bt = [b, None, b_fact[1]]
    
    m = b
    lt = 1
    t = b_fact[0]
    mr = MontgomeryReduction(m, b_fact)
    
    n = 0
    ps = [False, None]
    
    while x < max_pow and y < max_pow and not ps[0]:       
        
        x = get_x(a, x, c, z, mr, t, lt)
        y = get_y(a, x, c, z, b, lrab, max_yc)
        
        lt = t
        m *= b
        t *= b
        mr.set_modulus(m, bt)
        
        ps = prob_solution(a, x, b, y, c, z)
        n += 1
        
        #print("{!s}: ({!s}, {!s})".format(n, x, y))
        
    return [x, y, ps]

def get_x(a, x, c, z, mr, t, lt):
    
    ra, rc = mr.convert(a), mr.convert(c)
    rcz = mr.exp(rc, z)
    
    f = mr.exp(ra, t - x)
    g = mr.exp(ra, lt)
    h = mr.mult(f, rcz)
    
    p = mr.dlog(g, h)   
    
    x += p * lt
    
    return x

#find y by chinese remainder theorem modulo prime squares (allows chinese remainder theorem modulo primorial larger than y
def get_y(a, x, c, z, b, lrab, max_yc):
    
    y_res_mod = []
    
    max_ya = x * lrab
    max_y = max(max_ya, max_yc)
    
    p = 1
    mr = MontgomeryReduction()
    
    for i, pprime in enumerate(prob_primes):
        if is_gen_n_k(b, pprime, all_fact[pprime]): 
            m = pprime**2
            
            #squaring ensures totient is divisible by prime
            mr.set_modulus(m, prime_sq_fact[i])
                    
            ra, rc, rb = mr.convert(a), mr.convert(c), mr.convert(b)
            rax, rcz = mr.exp(ra, x), mr.exp(rc, z)
            
            rd = (rax - rcz) % m
            
            if is_coprime(rd, m):                
                y_sp = mr.dlog_mod(rb, rd, [pprime, 1])
                y_res_mod.append([y_sp, pprime])
                
                p *= pprime
                
                if p > max_y:
                    break
    else:
        ipdb.set_trace()
        print("Not enough primes to get y")
        #raise Error("Not enough primes to get y")
    
    return crt(y_res_mod)

def prob_solution(a, x, b, y, c, z):
    
    num_pp = len(safe_primes)
        
    mr = MontgomeryReduction()
    plus = True
    minus = True
    
    for i in range(num_ps_tests):
        #using safe primes ensures that the order of a and b modulo m is high, thus ensuring the probability of a^x+/-b^y=c^z is low, on the order of 1/m
        m = safe_primes[randrange(num_pp)]
        
        mr.set_modulus(m)
        
        ra, rb, rc = mr.convert(a), mr.convert(b), mr.convert(c)
        rax, rby, rcz = mr.exp(ra, x), mr.exp(rb, y), mr.exp(rc, z)
        
        if plus and (rax + rby) % m != rcz:
            plus = False
            s = '-'
            
            if not minus:
                break
        
        if minus and (rax - rby) % m != rcz:
            minus = False
            s = '+'
            
            if not plus:
                break
    else:
        
        return [True, s]
    
    return [False, s]

def populate_all(n=int_ceiling):
    
    global prob_primes, safe_primes, all_fact, prime_sq_fact
    
    prob_primes, safe_primes, all_fact, prime_sq_fact = get_all_lists(n)
