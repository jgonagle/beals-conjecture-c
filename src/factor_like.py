import os
import sys
from math import ceil
from time import strftime
from functools import partial
from multiprocessing import Lock, Process, Manager, Value, BoundedSemaphore, cpu_count
from sympy import symbols, sympify, poly
from memory_profiler import profile
from ipdb import set_trace

import populate_num_lists as pnl
from time_wrap import timed
from num_theory import is_mutually_coprime

int_ceiling = 1000000
threshhold_gap = 100
serial_gap = 1

prob_primes = []
semiprimes = []
primes_semiprimes = []
all_fact = []
psp_fact = {}

log_dir = "../data/logs/"
sol_dir = "../data/solutions/"

log_path = None
sol_path = None

print_lock = None
log_lock = None
sol_lock = None

a, b = symbols("a b")

@timed
#@profile
def find_all_ab(psp, min_xy=3, max_xy=10, min_c=3, max_c=int_ceiling, min_z=3, max_z=10):

    set_log_path(psp, min_xy, max_xy, min_c, max_c, min_z, max_z)
    set_sol_path(psp, min_xy, max_xy, min_c, max_c, min_z, max_z)
    
    set_new_print_lock()
    set_new_log_lock()
    set_new_sol_lock()
    
    if psp:
        populate_all_psp_lists(max_c)
    else:
        populate_all_lists(max_c)
    
    xyz_gen = spiral_sum(min_xy, max_xy, min_z, max_z, all=True)
    
    print_write(True, "Starting search for xy in [{!s},{!s}), {} c in [{!s},{!s}), and z in [{!s},{!s})".format(min_xy, max_xy, "psp" if psp else "all", min_c, max_c, min_z, max_z))
    
    for xy, z in xyz_gen:
        if xy != z:
            search_c_parallel(psp, xy, xy + 1, min_c, max_c, z, z + 1)
            
    print_write(True, "\nCompleted search for xy in [{!s},{!s}), {} c in [{!s},{!s}), and z in [{!s},{!s})".format(min_xy, max_xy, "psp" if psp else "all", min_c, max_c, min_z, max_z))
    
#@timed
#@profile
def search_c_parallel(psp, min_xy=3, max_xy=10, min_c=3, max_c=int_ceiling, min_z=3, max_z=10):
    
    serial_jobs = get_job_gen(min_xy, max_xy, min_c, max_c, min_z, max_z)
    
    sem_size = cpu_count()
    proc_sem = BoundedSemaphore(sem_size)
    max_tested = Value('i', -1)
    next_threshhold = Value('i', threshhold_gap)

    def sync_search(search_range):
        
        search_c_serial(psp, *search_range)
        proc_sem.release()
        
        with max_tested.get_lock():
            max_tested.value = max(max_tested.value, search_range[3])
            
            #unnecessary; max_tested lock takes care of atomicity. included for maintainability
            with next_threshhold.get_lock():
                if max_tested.value > next_threshhold.value:
                    print_write(True, "\t\t{!s}".format(next_threshhold.value))
                    next_threshhold.value += threshhold_gap
    
    print_write(True, "\n\tStarting search for xy in [{!s},{!s}), {} c in [{!s},{!s}), and z in [{!s},{!s})\n".format(min_xy, max_xy, "psp" if psp else "all", min_c, max_c, min_z, max_z))
    
    for s_job in serial_jobs:
        proc_sem.acquire()        
        p = Process(target=sync_search, args=(s_job,))
        p.start()
        
    wait_sem(proc_sem, sem_size)
    
    print_write(True, "\n\tCompleted search for xy in [{!s},{!s}), {} c in [{!s},{!s}), and z in [{!s},{!s})".format(min_xy, max_xy, "psp" if psp else "all", min_c, max_c, min_z, max_z))

#@timed
#@profile
def search_c_serial(psp, min_xy=3, max_xy=10, min_c=3, max_c=int_ceiling, min_z=3, max_z=10):
    
    #print_write(True, "\n\t\tStarting search for xy in [{!s},{!s}), {} c in [{!s},{!s}), and z in [{!s},{!s})".format(min_xy, max_xy, "psp" if psp else "all", min_c, max_c, min_z, max_z))
    
    #print_write(True, "\t\t{!s}-{!s}".format(min_c, max_c))
    
    if psp:
        cl = filter(lambda x: min_c <= x < max_c, primes_semiprimes)
    else:
        cl = range(min_c, max_c)
    
    for xy in range(min_xy, max_xy):
        pp = get_p(xy, True)
        pm = get_p(xy, False)
        
        for z in range(min_z, max_z):
            
            for c in cl:
                if psp:
                    czd = get_all_psp_divisors(c, z)
                else:
                    czd = get_all_comp_divisors(c, z)
                
                cz = czd[-1]
                
                for clp in czd:
                    #techincally, we could just use its mirror image since the divisor list is sorted and a[i] * a[len(a) - (i + 1)] = c^z
                    crp = cz // clp
                    sm = solve(pm, clp, crp, False)
                    
                    for a in sm:
                        b = a - clp
                        
                        if b > 0 and is_mutually_coprime(a, b, c) and b**xy + cz == a**xy:
                            print_write(False, "{!s}^{!s} + {!s}^{!s} = {!s}^{!s}".format(b, xy, c, z, a, xy))
                            
                    if xy % 2 == 1:
                        sp = solve(pp, clp, crp, True)
                        
                        for a in sp:
                            b = clp - a
                            
                            if b > 0 and is_mutually_coprime(a, b, c) and a**xy + b**xy == cz:
                                print_write(False, "{!s}^{!s} + {!s}^{!s} = {!s}^{!s}".format(a, xy, b, xy, c, z))
    
    #print_write(True, "\n\t\tCompleted search for xy in [{!s},{!s}), {} c in [{!s},{!s}), and z in [{!s},{!s})".format(min_xy, max_xy, "psp" if psp else "all", min_c, max_c, min_z, max_z))
                    
def get_job_gen(min_xy, max_xy, min_c, max_c, min_z, max_z):
         
    for c in range(min_c, max_c, serial_gap):
        yield [min_xy, max_xy, c, min(max_c, c + serial_gap), min_z, max_z]

#@timed
#@profile
def solve(p, l, r, plus):
    
    q = get_q(l, plus)
    f = p.subs(b, q)
    v = poly(f - r).intervals()
    
    s = []
    
    if v:
        s = [i[0][0] for i in v if i[0][0] > 0 and i[0][0] == i[0][1]]
    
    return s

def get_p(n, plus):
    
    psa = ['' for t in range(n)]
    rn = n - 1

    for i in range(0, rn + 1):
        if i == 0:
            s = ''
        else:
            if plus:
                if i % 2 == 0:
                    s = '+'
                else:
                    s = '-'
            else:
                s = '+'

        psa[i] = ''.join([s, "a**", str(rn - i), '*', "b**", str(i)])

    ps = ''.join(psa)
    p = sympify(ps)
    
    return p

def get_q(d, plus):
    
    if plus:
        qs = ''.join([str(d), '-a'])
    else:
        qs = ''.join(['a-', str(d)])
    
    q = sympify(qs)
    
    return q

def get_all_comp_divisors(n, re):
    
    fact = all_fact[n][1]
    fact_e = [[p, re * e] for p, e in fact]
    
    ad = divisor_recursive(fact_e)

    return sorted(ad)

def get_all_psp_divisors(psp, re):
    
    fact = psp_fact[psp]
    fact_e = [[p, re * e] for p, e in fact]
    
    ad = divisor_recursive(fact_e)

    return sorted(ad)

def divisor_recursive(fact_e):
    
    if not fact_e:
        return [1]
    else:
        p, e = fact_e[0]
        dr = divisor_recursive(fact_e[1:])
        
        ad = []
        m = 1
        
        for i in range(e + 1):
            ad.extend([m * d for d in dr])
            
            m *= p
        
        return set(ad) 

#@timed
def populate_all_psp_lists(n=int_ceiling):
    
    global primes_semiprimes, psp_fact
    
    primes_semiprimes, psp_fact = pnl.get_all_psp_lists(n)

#@timed
def populate_all_lists(n=int_ceiling):
    
    global prob_primes, semiprimes, primes_semiprimes, all_fact, psp_fact
    
    prob_primes, semiprimes, primes_semiprimes, safe_primes, all_fact, psp_fact, prime_sq_fact = pnl.get_all_lists(n)

def spiral_sum(min_x=0, max_x=int_ceiling, min_y=0, max_y=int_ceiling, n=int_ceiling, all=False):
    
    c = 0
    
    if min_x <= 0 < max_x and min_y <= 0 < max_y:
        yield [0, 0]
        c += 1
    
    #start origin
    x = 0
    y = 1
    
    xi = 1
    yi = -1
    
    while all or c < n:    
        if min_x <= x < max_x and min_y <= y < max_y:
            yield [x, y]
            c += 1
        
        if x == max_x - 1 and y == max_y - 1:
            break
        
        if x + xi == 0 and y + yi > 0:     
            x = 0
            y += 2
            
            yi = -1
        else:
            x += xi
            y += yi
            
            if x > 0 and y == 0:
                xi = -1
            elif x == 0 and y < 0:
                yi = 1
            elif x < 0 and y == 0:
                xi = 1

def print_write(is_log, some_str):
    
    if is_log:
        file_lock = log_lock
        file_path = log_path
    else:
        file_lock = sol_lock
        file_path = sol_path        
        
    file_lock.acquire()
    with open(file_path, 'a') as file_writer:
        file_writer.write(some_str + '\n')
    file_lock.release()
    
    if not is_log:
        some_str = "\t\t\t" + some_str
    
    print_lock.acquire()
    print(some_str)
    print_lock.release()
    

def set_log_path(psp, min_xy, max_xy, min_c, max_c, min_z, max_z):
    
    global log_path
    
    filename = "{}{}_xy({!s},{!s})_c({!s},{!s})_z({!s},{!s})_{}".format(log_dir, "psp" if psp else "all", min_xy, max_xy, min_c, max_c, min_z, max_z, get_time_str())
    
    local_path = os.path.abspath(os.path.dirname(__file__))
    log_path = os.path.join(local_path, filename)
    
def set_sol_path(psp, min_xy, max_xy, min_c, max_c, min_z, max_z):
    
    global sol_path
    
    filename = "{}{}_xy({!s},{!s})_c({!s},{!s})_z({!s},{!s})_{}".format(sol_dir, "psp" if psp else "all", min_xy, max_xy, min_c, max_c, min_z, max_z, get_time_str())
    
    local_path = os.path.abspath(os.path.dirname(__file__))
    sol_path = os.path.join(local_path, filename)
    
def set_new_print_lock():
    
    global print_lock
    
    print_lock = Lock()

def set_new_log_lock():
    
    global log_lock
    
    log_lock = Lock()

def set_new_sol_lock():
    
    global sol_lock
    
    sol_lock = Lock()
    
def wait_sem(sem, sem_size):
    """blocks until semaphore count is zero
    
    acquires number of semaphore slots equal to the semaphore's size and then releases all slots.
    sem must be a sempahore and count must be equal to the semaphore's size.  a count greater than the semaphore's size will result in an error."""
    
    
    for i in range(sem_size):
        sem.acquire()
    
    for i in range(sem_size):
        sem.release()

def get_time_str():
    
    return strftime("%b_%d_%Y_%H:%M:%S")

if __name__ == "__main__":
    arguments = sys.argv[1:]
    
    psp = True if arguments[0] == "True" else False
    min_xy = int(arguments[1])
    max_xy = int(arguments[2])
    min_c = int(arguments[3])
    max_c = int(arguments[4])
    min_z = int(arguments[5])
    max_z = int(arguments[6])
    
    find_all_ab(psp, min_xy, max_xy, min_c, max_c, min_z, max_z)