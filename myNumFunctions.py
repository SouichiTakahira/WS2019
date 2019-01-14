import numpy as np
import random

# -----------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------

def get_random_unit_vec(dim):
    vec = []
    for i in range(dim):
        vec.append(complex(random.uniform(-1,1), random.uniform(-1,1)))
    return vec/(np.linalg.norm(vec))

def normalize(vec):
    return vec/(np.linalg.norm(vec))

def get_standard_basis(dim, k):
    vec = []
    for i in range(dim):
        vec.append(complex(0,0))
    vec[k] = complex(1,0)
    return vec

def Rx(t):
    return np.matrix([[np.cos(t/2),-np.sin(t/2)*1j],[-np.sin(t/2)*1j,np.cos(t/2)]])

def Ry(t):
    return np.matrix([[np.cos(t/2),-np.sin(t/2)],[np.sin(t/2),np.cos(t/2)]])

def Rz(t):
    return np.matrix([[np.exp(-t*1j/2), 0],[0, np.exp(t*1j/2)]])

def Id(t):
    return np.matrix([[np.exp(-(t/2)*1j),0],[0,np.exp(-(t/2)*1j)]])
