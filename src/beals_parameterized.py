import math
import sys

def is_root(r, s, p):
    
    c = (64 * r**9 * s**3 + 3 * r**8 + 576 * r**7 * s**5 - 36 * r**6 * s**2 +
          1728 * r**5 * s**7 + 162 * r**4 * s**4 + 1728 * r**3 * s**9 -
          324 * r**2 * s**6 + 243 * s**8)
    d = 1
    
    if c > 0 and c % d == 0:
    
        c = c // d
        
        e = get_root(c, p)
    
        return [r**e == c, e]

    else:
        return [False, None]
    
def chakravala(n):
    
    m = 1
    k = 1
    
    a = 1
    b = 0
 
    while k != 1 or b == 0:
        m = k * ( m // k + 1) - m
        m = m - ((m - get_root(n, 2) - 1) // k) * k
 
        t = (a * m + n * b) // abs(k)
        b = (a + b * m) // abs(k)
        k = (m**2 - n) // k
 
        a = t
 
    return [a, b] 

def get_root(n, p):
    
    s, e = find_search_range(n, p)
    
    return root_binary_search(s, e, n, p)

def find_search_range(n, p):

    l = 1
    
    #establish limits of search for binary search for root(n)
    while (l**p < n):
        l *= 2
        
    s = l // 2
    e = l

    return [s, e]

#returns smallest integer whose square is greater than or equal to n
def root_binary_search(s, e, n, p):
    
    while (s <= e):
        m = (s + e) // 2
        
        if (m**p < n):
            s = m + 1
        elif (m**p > n):
            e = m - 1
        else:
            return m
    
    return s - 1

def get_all_nth_pow_solutions(all=False, nn=50):
    
    n = 0
    e = 2
    
    a = 0
    b = 2
    
    i = 1
    
    while(all or n < nn):
        
        try:
            n += 1
        
            a += i
            b += -i
            
            if a == 0 or b == 0:
                
                a -= min(0, i)
                b -= min(0, -i)
                
                i *= -1
            
            for j in [-1, 1]:
                for k in [-1, 1]:
            
                    g = is_root(j * a, k * b, e)
                    #print('\t', j * a, ',', k * b, ',', g[1])

                    if g[0]:
                        print("!!!!", j * a, ',', k * b, ',', g[1], "!!!!")
                        
                        raise KeyboardInterrupt
        
        except KeyboardInterrupt: 
            print("\nn:", n)
            print("a, b, c:", a, ',', b, ',', g[1])
            
            yes = set(["yes","ye", 'y', ''])
            no = set(["no",'n'])
            
            choice = '?'
            
            #until a valid choice has been made, prompt for a yes or no answer on whether to continue
            while (choice not in yes) and (choice not in no):
            
                choice = input("\nContinue? (Y/N):").lower()
                
            if choice in no:
                sys.exit(0)
                
def find_square_chak(all=False,nn=50):
    
    n = 0
    e = 2
    
    while(all or n < nn):
        try:                
            n += 1
            
            c = chakravala(3 * n**4)[1]
            g = is_root(c, n, e)

            if g[0]:
                print("!!!!", c, ':', g[1], "!!!!")
                
                raise KeyboardInterrupt
        
        except KeyboardInterrupt:            
            print("\nn:", n,)
            print(c, ',', g[1])
               
            yes = set(["yes","ye", 'y', ''])
            no = set(["no",'n'])
            
            choice = '?'
            
            #until a valid choice has been made, prompt for a yes or no answer on whether to continue
            while (choice not in yes) and (choice not in no):
               
                choice = input("\nContinue? (Y/N):").lower()
                
            if choice in no:
                sys.exit(0)

def get_all_chak_powers(chak, d, all=False, nn=50):
    
    n = 0
    e = 2
    
    x, y = chakravala(chak)
    
    cx , cy = x, y
    
    g = is_root(cx, cy, d, e)

    if g[0]:
        print("!!!!", cy, ':', g[1], "!!!!")
        
    g = is_root(-cx, cy, d, e)

    if g[0]:
        print("!!!!", cy, ':', g[1], "!!!!")
    
        raise Keyboard_Interrupt
    
    while(all or n < nn):
        try:                
            n += 1

            tx = cx * x + chak * cy * y
            ty = cx * y + cy * x
            
            cx = tx
            cy = ty
            
            g = is_root(cx, cy, d, e)
            
            #print(cx, cy)

            if g[0]:
                print("!!!!", cy, ':', g[1], "!!!!")
                
                raise KeyboardInterrupt
            
            g = is_root(-cx, cy, d, e)

            if g[0]:
                print("!!!!", cy, ':', g[1], "!!!!")
        
        except KeyboardInterrupt:            
            print("\nn:", n,)
            print(cx, ',',       cy)
               
            yes = set(["yes","ye", 'y', ''])
            no = set(["no",'n'])
            
            choice = '?'
            
            #until a valid choice has been made, prompt for a yes or no answer on whether to continue
            while (choice not in yes) and (choice not in no):
               
                choice = input("\nContinue? (Y/N):").lower()
                
            if choice in no:
                sys.exit(0)
                
def get_all_beals_solutions(t, n):
    
    for i in range(3, t, 2):
        get_all_pow_solutions(i, nn=n, all=False)
        
def zig_zag(all=False, nn=50):
    
    n = 0
    
    a = 1
    b = 1
    
    c = 1
    
    print(a, ',', b)
    
    while(all or n < nn):
        
        n += 1
        
        a += c
        b += -c
        
        if a == 0 or b == 0:
            
            a -= min(0, c)
            b -= min(0, -c)
            
            c *= -1
        
        print(a, ',', b)
            