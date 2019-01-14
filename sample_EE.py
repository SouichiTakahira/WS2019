import matplotlib.pyplot as plt
from qiskit import ClassicalRegister, QuantumRegister
from qiskit import QuantumCircuit, execute, BasicAer

from numpy import *
from myNumFunctions import *
from myPrintFunctions import *
import random
import cmath
import sys

# ----------------------------------------------------------------------
# Hadamard gates

def hadamards(circ, qreg):
    for j in range(qreg.size):
        circ.h(qreg[j])

# ----------------------------------------------------------------------
# Controlled expiAj operation

def ctrl_expiAtpow(circ, ctrl_qreg, targ_qreg, A):
    for j in range(ctrl_qreg.size):
        for k in range(2**(ctrl_qreg.size-j-1)):
        	ctrl_expiAt(circ, ctrl_qreg[j], targ_qreg[0], A, -pi)

# ----------------------------------------------------------------------
# The Hamiltonian simulation algorithm  
# 
# Nielsen, Michael A., and Isaac L. Chuang. 
# "Quantum computation and quantum information." (2000).
# Exercise 4.50. 

def ctrl_expiAt(qc, ctrl, targ, A, dt):
    a = A[0,0].real
    b = A[1,1].real
    c = A[1,0].real
    d = A[1,0].imag
    theta_a = a*dt
    theta_b = b*dt
    theta_c = c*dt
    theta_d = d*dt
    ctrl_phase(theta_a, ctrl, targ, qc)
    qc.crz(theta_a, ctrl, targ)
    ctrl_phase(theta_b, ctrl, targ, qc)
    qc.crz( -theta_b, ctrl, targ)
    crx(2*theta_c, ctrl, targ, qc)
    cry(2*theta_d, ctrl, targ, qc)
    cry(2*theta_d, ctrl, targ, qc)
    crx(2*theta_c, ctrl, targ, qc)
    qc.crz( -theta_b, ctrl, targ)
    ctrl_phase(theta_b, ctrl, targ, qc)
    qc.crz(  theta_a, ctrl, targ)
    ctrl_phase(theta_a , ctrl, targ, qc)

def ctrl_expiAt_INV(qc, ctrl, targ, A, dt):
    ctrl_expiAt(qc, ctrl, targ, A, -dt)

def ctrl_phase(lam, ctrl, targ, qc):
    qc.cu1(-lam/2, ctrl, targ)
    qc.cx(ctrl, targ)
    qc.cu1(-lam/2, ctrl, targ)
    qc.cx(ctrl, targ)

def crx(lam, ctrl, targ, qc):
    qc.cu3(lam, -pi/2, pi/2, ctrl, targ) # u3(lam, -pi/2, pi/2) = Rx(lam)

def cry(lam, ctrl, targ, qc):
    qc.cu3(lam, 0, 0, ctrl, targ) # u3(lam, 0, 0) = Ry(lam)

# ----------------------------------------------------------------------
# The inverse quantum Fourier transform 

def qft_INV(circ, qreg):
    for j in range(qreg.size):
        circ.h(qreg[j])
        for k in range(2, qreg.size-j+1):
            circ.cu1(-2*pi/(2**k), qreg[j+k-1], qreg[j])
    reverse(circ, qreg)

def reverse(circ, qreg):
    for j in range(qreg.size // 2):
        circ.swap(qreg[j], qreg[qreg.size-1-j])

# ----------------------------------------------------------------------
# eigenvalues estimation

def eigenvalues_estimation(circ, ctrl, targ, A):
    hadamards(circ, ctrl)
    ctrl_expiAtpow(circ, ctrl, targ, A)
    qft_INV(circ, ctrl)
# ----------------------------------------------------------------------
# main

def main():
    A = array([[3/8,1/4],[1/4,3/8]])
    psi_vec = array([1.0, 0.0])
    myinfo(A, psi_vec)
   
    m = int(sys.argv[1])
    qr = QuantumRegister(m, 'qr')
    cr = ClassicalRegister(m, 'cr')
    psi = QuantumRegister(1, 'psi')
    qc = QuantumCircuit(qr, psi, cr)

    qc.initialize(psi_vec, psi) 
    eigenvalues_estimation(qc, qr, psi, A)

    for j in range(m):
        qc.measure(qr[j], cr[j])

    backend = BasicAer.get_backend('qasm_simulator')
    job = execute(qc, backend, shots=1024)
    result = job.result()
    counts = result.get_counts(qc)
    print(counts)  

# ----------------------------------------------------------------------
# printout the information of the matrix A and the vector |psi>

def myinfo(M, psi):
    eigs,eigenvecs = linalg.eig(M)
    print("------------------------------------------")
    print("INFORMATION OF THE SETTING")
    print("------------------------------------------")
    print("coefficient A:\n", M)
    print("------------------------------------------")
    print("0th eigenvector |u0> :", eigenvecs[:,0])
    print("1th eigenvector |u1> :", eigenvecs[:,1])
    print("------------------------------------------")
    print("0th eigenvalue  lam0 :", eigs[0].real)
    print("1th eigenvalue  lam1 :", eigs[1].real)
    print("------------------------------------------")
    print("|psi>:", psi)

main()
