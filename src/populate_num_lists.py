import math
from random import randrange
import ipdb
from time_wrap import timed
from num_theory import exp_bs

int_ceiling = 10000
prob_comp_ind = .001
prob_comp_list = .1

prob_primes = []
semiprimes = []
primes_semiprimes = []
safe_primes = []
all_fact = []
psp_fact = {}
prime_sq_fact = []

def get_all_lists(n=int_ceiling):
    
    populate_all_lists(n)
    
    return [prob_primes, semiprimes, primes_semiprimes, safe_primes, all_fact, psp_fact, prime_sq_fact]

def get_all_psp_lists(n=int_ceiling):
    
    populate_all_psp_lists(n)
    
    return [primes_semiprimes, psp_fact]

#@timed
def get_pprimes_list(n=int_ceiling):
    
    populate_prob_primes(n)
    
    return prob_primes

def populate_all_lists(n=int_ceiling):

    populate_prob_primes(n)
    populate_semiprimes(n, repopulate_others=False)
    populate_primes_semiprimes(n, repopulate_others=False)
    populate_safe_primes(n, repopulate_others=False)
    populate_all_fact(n, repopulate_others=False)
    populate_psp_fact(n, repopulate_others=False)
    populate_prime_sq_fact(n, repopulate_others=False)
    
def populate_all_psp_lists(n=int_ceiling):
    
    populate_prob_primes(n)
    populate_semiprimes(n, repopulate_others=False)
    populate_primes_semiprimes(n, repopulate_others=False)
    populate_psp_fact(n, repopulate_others=False)

def populate_prob_primes(n=int_ceiling):
    
    global prob_primes
    
    #miller-rabin test doesn't work on 2 or 3 (empty random int range)
    prob_primes = [2, 3]
    
    #want probability of there being a probable prime in our list to be < prob_comp_list
    pci = 4**(math.log(prob_comp_list, 4) + math.log(1 / n, 4))
        
    for i in range(5, n, 2):
        if miller_rabin_test(i, pci):
            prob_primes.append(i)
            
def populate_semiprimes(n=int_ceiling, repopulate_others=True):
    """repopulates the semiprimes list by getting all semiprimes up to n
    
    repopulate_others only to be False if we're sure prob_primes contains all primes up to n.  Otherwise, semiprimes will contain semiprimes for all primes up to the prob_primes max"""
    
    global semiprimes, prob_primes
    
    semiprimes = []
    
    #no side effects, other than safe_primes
    tpp = prob_primes
    
    if repopulate_others:
        populate_prob_primes(n)
    
    for i, prime_a in enumerate(prob_primes):
        for prime_b in prob_primes[i + 1:]:
            sp = prime_a * prime_b
            
            if sp >= n:
                break
            
            semiprimes.append(sp)
                
    semiprimes = sorted(semiprimes)
     
     #no side effects, other than semiprimes
    prob_primes = tpp
 
def populate_primes_semiprimes(n=int_ceiling, repopulate_others=True):
    """repopulates the primes_semiprimes list by getting all primes and semiprimes up to n
    
    repopulate_others only to be False if we're sure prob_primes and semiprimes contains all primes and semiprimes, respectively, up to n.  Otherwise, primes_semiprimes will contain primes and semiprimes up to the prob_primes max and semiprimes max, respectively"""
    
    global primes_semiprimes, prob_primes, semiprimes
        
    #no side effects, other than prime_semiprimes
    tpp = prob_primes
    tsp = semiprimes
    
    if repopulate_others:
        populate_prob_primes(n)
        populate_semiprimes(n, repopulate_others=False)
    
    primes_semiprimes = sorted(prob_primes + semiprimes)
    
    #no side effects, other than primes_semiprimes
    prob_primes = tpp
    semiprimes = tsp    

def populate_safe_primes(n=int_ceiling, repopulate_others=True):
    """repopulates the safe_primes list by getting all safe primes up to n
    
    repopulate_others only to be False if we're sure prob_primes contains all primes up to n.  Otherwise, safe_primes will contain safe_primes for all primes up to the prob_primes max"""
    
    global safe_primes, prob_primes
    
    safe_primes = []
    
    #no side effects, other than safe_primes
    tpp = prob_primes
    
    if repopulate_others:
        populate_prob_primes(n)
    
    pp_set = set(prob_primes)
    
    for pp in prob_primes:        
        if (pp - 1) // 2 in pp_set:
            safe_primes.append(pp)
    
    #no side effects, other than safe_primes
    prob_primes = tpp

#@timed            
def populate_all_fact(n=int_ceiling, repopulate_others=True):
    """repopulates the all_fact list by factoring all numbers up to n
    
    repopulate_others only to be False if we're sure prob_primes contains all primes up to n.  Otherwise, we'll get caught in an infinite loop below.  Could be alleviated by breaking from for loop over primes instead of using while loop.  Even still, all numbers less than n will still not be fully factored"""
    
    global all_fact, prob_primes
    
    #no side effects, other than all_fact
    tpp = prob_primes
    
    if repopulate_others:
        populate_prob_primes(n)
    
    pprimes_set = set(prob_primes)
    
    all_fact = [None for i in range(n)]
    
    #techincally don't need, but it makes accessing an number's factorization and totient simpler
    #format is totient, factorization, and totient factorization
    all_fact[0] = [0, [[0, 1]], [[0, 1]]]
    all_fact[1] = [0, [[1, 1]], [[0, 1]]]
    
    for i in range(2, n):
        
        if i in pprimes_set:
            all_fact[i] = [i - 1, [[i, 1]], all_fact[i - 1][1]]
        else:        
            d = i
            p_n = 0
            
            tot = 1
            fact = []
            tot_fact = None
            
            while d != 1:
                p = prob_primes[p_n]
                c = 0
                
                while d % p == 0:
                    c += 1
                    d //= p
                
                if c > 0:
                    tot *= (p - 1) * p**(c - 1)
                    fact.append([p, c])
                
                p_n += 1
            
            #don't need to calculate the totient factorization since totient < i
            all_fact[i] = [tot, fact, all_fact[tot][1]]
    
    #no side effects, other than all_fact
    prob_primes = tpp
 
#@timed
def populate_psp_fact(n=int_ceiling, repopulate_others=True):
    """repopulates the psp_fact list by getting the factorizations of all primes and semiprimes up to n
    
    repopulate_others only to be False if we're sure prob_primes and primes_semiprimes contains all primes and semiprimes, respectively, up to n.  Otherwise, p"_sp_fact will contain primes and semiprimes up to the prob_primes max and primes_semiprimes max"""
    
    global psp_fact, prob_primes, primes_semiprimes
    
    #no side effects, other than psp_fact
    tpp = prob_primes
    tpsp = primes_semiprimes
    
    if repopulate_others:
        populate_prob_primes(n)
        populate_semiprimes(n, repopulate_others=False)
        populate_primes_semiprimes(n, repopulate_others=False)
    
    psp_fact = {}
    pprimes_set = set(prob_primes)    
    
    for psp in primes_semiprimes:
        if psp in pprimes_set:
            psp_fact[psp] = [[psp, 1]]
        else:
            for p in prob_primes:
                if psp % p == 0:
                    psp_fact[psp] = [[p, 1], [psp // p, 1]]
                    
                    break
                
    #no side effects, other than psp_fact
    prob_primes = tpp
    primes_semiprimes = tpsp
                        
def populate_prime_sq_fact(n=int_ceiling, repopulate_others=True):
    """repopulates the prime_sq_fact list by factoring all prime squares for all primes up to n
    
    repopulate_others only to be False if we're sure prob_primes contains all primes up to n.  Otherwise, prime_sq_fact will contain prime square factorizations for all primes up to the prob_primes max"""
    
    global prime_sq_fact, prob_primes, all_fact
    
    #no side effects, other than prime_sq_fact
    tpp = prob_primes
    taf = all_fact
    
    if repopulate_others:
        populate_prob_primes(n)
        populate_all_fact(n)
    
    prime_sq_fact = [None for i in range(len(prob_primes))]
    
    for i, pprime in enumerate(prob_primes):
        p_fact = all_fact[pprime]        
        prime_sq_fact[i] = [pprime * p_fact[0], [[pprime, 2]],  p_fact[2] + [[pprime, 1]]]
    
    #no side effects, other than prime_sq_fact
    prob_primes = tpp
    all_fact = taf
             
def miller_rabin_test(n, pci=prob_comp_ind):
    
    nm = n - 1
    d = nm
    s = 0
    
    #number of tests required to establish probable primeness    
    k = math.ceil(-math.log(pci, 4))
    
    while d % 2 == 0:
        s += 1
        d //= 2
        
    sm = s - 1
    
    for i in range(k):
        a = randrange(2, nm)
        
        x = exp_bs(a, d, n)
        
        if x == 1 or x == nm:
            continue
        else:
            for j in range(sm):
            
                x = x**2 % n
                
                #early stopping x=1 for all powers of s > j
                if x == 1:
                    return False
                
                #early stopping, the condition we're trying to prevent
                if x == nm:
                    break
            else:
                return False
    
    return True
