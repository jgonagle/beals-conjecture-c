import populate_num_lists as pnl
from time_wrap import timed
from num_theory import is_coprime, exp_bs

pprime_ceiling = 1000
prob_primes = []

def prob_beals(min_xyz=3, min_abc=2):
    
    populate_pprimes_list()
    
    #minimum sum is x+y+z=3+3+3=9 and a+b+c=2+3+5=10
    s = 19
    
    while True:
        #print("Shell with sum {!s} completed!".format(s))
        
        for a in range(min_abc, s + 1):
            s_1 = s - a
            
            for b in range(min_abc, s_1 + 1):
                s_2 = s_1 - b
                
                if b > a and is_coprime(a, b):         
                    
                    for c in range(min_abc, s_2 + 1):
                        s_3 = s_2 - c
                        
                        if is_coprime(a, c) and is_coprime(b, c):               
                            
                            for x in range(min_xyz, s_3 + 1):
                                s_4 = s_3 - x
                                
                                for y in range(min_xyz, s_4 + 1):
                                    z = s_4 - y
                                    
                                    if z >= min_xyz:
                                    
                                        for p in prob_primes:
                                            
                                            p_tot = p - 1
                                            
                                            a_m, b_m, c_m = a % p, b % p, c % p
                                            x_m, y_m, z_m = x % p_tot, y % p_tot, z % p_tot
                                            
                                            a_x, b_y, c_z = exp_bs(a, x, p), exp_bs(b, y, p), exp_bs(c, z, p)
                                            
                                            if (a_x + b_y) % p != c_z:
                                                break
                                        else:
                                            print("{!s}^{!s} + {!s}^{!s} = {!s}^{!s}".format(a, x, b, y, c, z))
                                        
        s += 1

#@timed
def populate_pprimes_list(n=pprime_ceiling):
    
    global prob_primes
    
    prob_primes = pnl.get_pprimes_list(n)