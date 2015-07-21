import math
import ipdb
from random import randrange
from functools import reduce
from time_wrap import timed
from populate_num_lists import get_all_lists
from num_theory import is_coprime, is_gen_n_k, odd_num_evens, crt, exp_bs, dlog, dlog_mod

int_ceiling = 10000
max_pow = 10**3

num_ps_tests = 5

prob_primes = []
safe_primes = []
all_fact = []
prime_sq_fact = []

def try_all_beals(lim=1000):
    
    populate_all()
    
    #b is in the primes since composite moduli don't have primitve roots.  Techincally, b could be in the prime powers or 2p^k, but the first is redundant for this application while the second isn't satisfied since b must be a power greater than two
    for b in prob_primes:
        
        print(b)        
        b_gp = get_gen_psq(b)
        
        for a in range(3, lim):
            if is_gen_n_k(a, b, all_fact[b]):            
                for c in range(3, lim):
                    if odd_num_evens(a, b, c) and is_coprime(a, c) and is_coprime(b, c):
                        for z in range(3, lim):
                            r = beals_solver(a, b, c, z, b_gp)
                            
                            #print(a, b, c, z)
                            
                            if r[2][0]:
                                if r[0] > 2 and r[1] > 2:
                                    print("\t!!!!{!s}^{!s} {} {!s}^{!s} = {!s}^{!s}!!!!".format(a, r[0], r[2][1], b, r[1], c, z))
                                    #wait = input("\tPossible solution found.  Press enter to continue.")
                                elif r[0] > 1 and r[1] > 1 and max(r[0], r[1]) > 2:
                                    print("\t{!s}^{!s} {} {!s}^{!s} = {!s}^{!s}".format(a, r[0], r[2][1], b, r[1], c, z))
                                    
#@timed
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
            min_y = tot_cover(b, num_af)
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
        
            b_gp = get_gen_psq(b)
            
            #print("\Solving {!s}^{!s} {!s} {!s}^{!s}...".format(a, x, s, b, y))
            r = beals_solver(a, b, c, 1, b_gp)
            
            assert r[0] == x and r[1] == y and r[2][0] and r[2][1] == s
            
            n_tests -= 1
            
    print("All tests passed!")            

@timed
def beals_solver(b, y, c, z, max_a, max_x):
    """finds a and x s.t. a^x + b^y = c^z
    
    b, c mutually coprime"""
    
    crt_mod = 1
    rep_crt = []
    
    for pprime in prob_primes:       
        
        pp_sq = pprime**2
        t = pprime * (pprime - 1)
        
        by = exp_bs(b, y, pp_sq)
        cz = exp_bs(c, z, pp_sq)
        d = (cz - by) % pp_sq
   
        res_exp_pairs = []
        
        for a in range(pp_sq):            
            l = 0 if a == 0 else 1
            
            for x in range(t):
                if l == d:
                    re_pair = [a % pprime, x % pprime]
                    
                    if re_pair not in res_exp_pairs:                        
                        res_exp_pairs.append(re_pair)
                else:
                    l = l * a % pp_sq
        
        rep_crt.append([pprime, res_exp_pairs])    
    
        crt_mod *= pprime
        
        if crt_mod > max_a and crt_mod > max_x:
            break
    
    np = reduce(lambda a, b: a * len(b[1]), rep_crt, 1)
    #np = reduce(lambda a, b: a * b[0], rep_crt, 1)
    print("{!s} chinese remainder theorem permutations".format(np))
    
    ax_pairs = crt_permutation(rep_crt, b, y, c, z)
    
    for axp in ax_pairs:
        
        a, x = axp
        
        print("\t{!s}^{!s} + {!s}^{!s} = {!s}^{!s}".format(a, x, b, y, c, z))

#TODO: prune leaves of crt permutation recursive tree so that a < max_a, x < max_x
    
def crt_permutation(rep_crt, b, y, c, z, crt_a=None, crt_x=None):
    
    if crt_a is None:
        crt_a = []
    
    if crt_x is None:
        crt_x = []
    
    if not rep_crt:
        a = crt(crt_a)
        x = crt(crt_x)
        
        ps = prob_solution(a, x, b, y, c, z)
        
        if ps[0]:
            return [[a, x]]
        else:
            return []
    
    car, cdr = rep_crt[0], rep_crt[1:]    
    ax_mod, res_exp_pairs = car[0], car[1]
    
    ax_pairs = []
    
    for re_pair in res_exp_pairs:
        a_res, x_res = re_pair        
        a_rm, x_rm = [[a_res, ax_mod]], [[x_res, ax_mod]]
        
        ax_pairs.extend(crt_permutation(cdr, b, y, c, z, crt_a + a_rm, crt_x + x_rm))
    
    return ax_pairs

def get_x(a, x, c, z, m, mf, t, lt):
    
    cz = exp_bs(c, z, m)
    
    f = exp_bs(a, t - x, m)
    g = exp_bs(a, lt, m)
    h = f * cz % m
    
    p = dlog(g, h, m, mf)
    
    x += p * lt
    
    return x

#find y by chinese remainder theorem modulo prime squares (allows chinese remainder theorem modulo primorial larger than y
def get_y(a, x, c, z, b, lrab, max_yc, b_gp):
    
    y_res_mod = []
    
    max_ya = x * lrab
    max_y = max(max_ya, max_yc)
    
    p = 1
    
    for pprime, pprime_sq, prime_sq_tot, prime_div in b_gp:
        #squared prime modulus ensures totient is divisible by prime
        ax, cz = exp_bs(a, x, pprime_sq), exp_bs(c, z, pprime_sq)
        
        d = (ax - cz) % pprime_sq
        
        if is_coprime(d, pprime_sq):                
            y_sp = dlog_mod(b, d, pprime_sq, prime_sq_tot, prime_div)
            y_res_mod.append([y_sp, pprime])
            
            p *= pprime
            
            if p > max_y:
                break
    else:
        print("Not enough primes to get y")
        #raise Error("Not enough primes to get y")
    
    return crt(y_res_mod)

def prob_solution(a, x, b, y, c, z):
    
    num_pp = len(safe_primes)
        
    plus = True
    minus = True
    
    for i in range(num_ps_tests):
        #using safe primes ensures that the order of a and b modulo m is high, thus ensuring the probability of a^x+/-b^y=c^z is low, on the order of 1/m
        m = safe_primes[randrange(num_pp)]
        
        ax, by, cz = exp_bs(a, x, m), exp_bs(b, y, m), exp_bs(c, z, m)
        
        if plus and (ax + by) % m != cz:
            plus = False
            s = '-'
            
            if not minus:
                break
        
        if minus and (ax - by) % m != cz:
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
    
def tot_cover(p, n):
    """calculates k s.t. phi(p^k) > n
    
    p must be prime"""
    return math.ceil(math.log(n / (p - 1), p)) + 1

def get_gen_psq(x):
    
    gen_prime_sq = []
    
    for i, pprime in enumerate(prob_primes):
        if is_gen_n_k(x, pprime, all_fact[pprime]):
            gen_prime_sq.append([pprime, pprime**2, prime_sq_fact[i][0], [pprime, 1]])

    return gen_prime_sq
