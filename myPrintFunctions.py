from numpy import *
import cmath

# ----------------------------------------------------------------------
# Printer

def print_state(state):
    dim = len(state)
    n   = int(log2(dim))
    for j in range(dim):
        print_amplitude_polar(j, state[j], n)
        
def print_cuted_state(state):
    dim = len(state)
    cut = 4
    n   = int(log2(dim))
    for j in range(0, cut):
        print_amplitude_polar(j, state[j], n)

def print_postselected_state(state):
    vec = array([state[1], state[3]])
    print_state(vec/linalg.norm(vec))

def print_amplitude_polar(index, amp, n):
    print(format(index), end='') # '0'+str(int(floor(log(dim)))-1))
    print(':', format(abs(amp), '0.10f'), end='')
    print(' *e^( '+format( cmath.phase(amp), ' 0.4f')+'i )', end='')
    print(' |'+format(index, '0'+str(n)+'b')+'> ')  

