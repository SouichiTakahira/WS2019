# Souichi Takahira
# Aichi Pref. Univ. 2018/1/14

import matplotlib.pyplot as plt
from qiskit import ClassicalRegister, QuantumRegister
from qiskit import QuantumCircuit, execute, BasicAer

from numpy import *
from myNumFunctions import *
from myPrintFunctions import *
import random
import cmath


# ----------------------------------------------------------------------
# Hadamard gates

def hadamards(circ, qreg):
    for j in range(qreg.size):
        circ.h(qreg[j])

# ----------------------------------------------------------------------
# Controlled U^j operation

def ctrl_iAt(circ, ctrl_qbit, targ_qreg, A):    
    ctrl_expiAt(circ, ctrl_qbit, targ_qreg[0], A, -pi)

def ctrl_iAt_INV(circ, ctrl_qbit, targ_qreg, A): 
    ctrl_expiAt(circ, ctrl_qbit, targ_qreg[0], A, pi)

# ---

def ctrl_iAtpow(circ, ctrl_qreg, targ_qreg, A):
    for j in range(ctrl_qreg.size):
        for k in range(2**(ctrl_qreg.size-j-1)):
            ctrl_iAt(circ, ctrl_qreg[j], targ_qreg, A)

def ctrl_iAtpow_INV(circ, ctrl_qreg, targ_qreg, A):
    for j in reversed(range(ctrl_qreg.size)):
        for k in reversed(range(2**(ctrl_qreg.size-j-1))):
            ctrl_iAt_INV(circ, ctrl_qreg[j], targ_qreg, A)

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
# Quantum Fourier transform and its inverse

def qft(circ, qreg):
    for j in range(qreg.size):
        circ.h(qreg[j])
        for k in range(2, qreg.size-j+1):
            circ.cu1(2*pi/(2**k), qreg[j+k-1], qreg[j])
    reverse(circ, qreg)

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
# Phase estimation algorithm

def eigenvalues_estimation(circ, ctrl, targ, A):
    hadamards(circ, ctrl)
    ctrl_iAtpow(circ, ctrl, targ, A)
    qft_INV(circ, ctrl)

def eigenvalues_estimation_INV(circ, ctrl, targ, A):
    qft(circ, ctrl)
    ctrl_iAtpow_INV(circ, ctrl, targ, A)
    hadamards(circ, ctrl)

# ----------------------------------------------------------------------
# Controlled rotations

def mult_cnot(circ, qreg, anci, mode): # mode: 0は白丸．1は黒丸
    ctrl_select(circ, qreg, mode)
    circ.ccx(qreg[0], qreg[1], anci[0])
    for j in range(1, qreg.size):
        circ.ccx(qreg[j], anci[j-1], anci[j])
    ctrl_select(circ, qreg, mode)

def mult_cnot_INV(circ, qreg, anci, mode):
    ctrl_select(circ, qreg, mode)
    for j in reversed(range(1, qreg.size)):
        circ.ccx(qreg[j], anci[j-1], anci[j])
    circ.ccx(qreg[0], qreg[1], anci[0])
    ctrl_select(circ, qreg, mode)

def ctrl_select(circ, qreg, mode):
    for j in range(len(mode)):
        if mode[j] == '0':
            circ.x(qreg[j])

def ctrl_rot_amp(circ, ctrl_qbit, targ_qbit, amp): ## amp <= 1
    theta = 2*arcsin(abs(amp))
    phase = cmath.phase(amp)
    circ.cu3(theta, 0, 0, ctrl_qbit, targ_qbit)
    circ.cu1(phase, ctrl_qbit, targ_qbit)

def ctrl_rot_one_routine(circ, ctrl_qreg, ctrl_anci, targ_qbit, amp, mode):
    mult_cnot(circ, ctrl_qreg, ctrl_anci, mode)
    ctrl_rot_amp(circ, ctrl_anci[ctrl_anci.size-1], targ_qbit, amp)
    mult_cnot_INV(circ, ctrl_qreg, ctrl_anci, mode)

def ctrl_rot(circ, ctrl_qreg, ctrl_anci, targ_qbit):
    n = ctrl_qreg.size
    for j in range(1, 2**n):
        mode = format(j, '0'+str(n)+'b')
        ctrl_rot_one_routine(circ, ctrl_qreg, ctrl_anci, targ_qbit, 1/j, mode)

# ----------------------------------------------------------------------
# HHL algorithm (original version)

def hhl_algorithm(circ, ctrl_qreg, ctrl_anci, anci_qbit, ket_b, A):
    eigenvalues_estimation(circ, ctrl_qreg, ket_b, A)
    ctrl_rot(circ, ctrl_qreg, ctrl_anci, anci_qbit)
    eigenvalues_estimation_INV(circ, ctrl_qreg, ket_b, A)


# ----------------------------------------------------------------------
# Measurement

def measurement(circ, qreg, creg):
    for j in range(qreg.size):
        circ.measure(qreg[j], creg[j])

# ----------------------------------------------------------------------
# Main 

def main():

    eps = 0
    A = array([[3/8+eps,1/4],[1/4,3/8+eps]])
    stateb = array( get_random_unit_vec(2)  )
    stateb = stateb/linalg.norm(stateb)
    myinfo(A, stateb)


    m = 4
    qr   = QuantumRegister(m, 'qr')
    qr_a = QuantumRegister(m, 'qr_a')
    anci = QuantumRegister(1, 'anci')
    ketb = QuantumRegister(1, 'ketb')

    cr   = ClassicalRegister(m, 'cr')
    cr_a = ClassicalRegister(1, 'cr_a')

    # |qr_a>|qr> o |b>|0>_a
    qc = QuantumCircuit(anci, ketb, qr, qr_a)

    qc.initialize(stateb, ketb)
    hhl_algorithm(qc, qr, qr_a, anci, ketb, A)

    get_result(qc, 0, m)

    # print(qc.draw(filename='test', output='text'))

def get_result(circ, mode, n):
    print("\nRESULT OF THE SIMULATOR:")
    if mode == 0:
        backend = BasicAer.get_backend('statevector_simulator')
        job = execute(circ, backend)
        result = job.result()
        statevector = result.get_statevector(circ)
        #
        # print(statevector[0:20])
        # myPrintFunction.pyにある関数を使用
        # print_state(statevector[0:20]) 
        print_postselected_state(statevector)
    else:
        backend = BasicAer.get_backend('qasm_simulator')
        job = execute(circ, backend, shots=128)
        result = job.result()
        counts = result.get_counts(circ)
        print(counts)  

    # print(circ)
    # print(circ.draw(scale=0.5, interactive=True, filename='testss', output='mpl'))

# ----------------------------------------------------------------------
# Information

def myinfo(M, b):
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
    print("|b>:", b)
    print("------------------------------------------")
    ip0 = vdot(eigenvecs[:,0], b)
    ip1 = vdot(eigenvecs[:,1], b)
    # print("<u0|b>:", ip0)
    # print("<u1|b>:", ip1)
    # print("------------------------------------------")
    vec0 = array([sqrt(1 - (1/(8*eigs[0].real))**2), 1/(8*eigs[0].real)])
    vec1 = array([sqrt(1 - (1/(8*eigs[1].real))**2), 1/(8*eigs[1].real)])
    vec = kron(ip0*eigenvecs[:,0], vec0) + kron(ip1*eigenvecs[:,1], vec1)
    print("<u0|b>(C/lam0)|u0> + <u1|b>(C/lam1)|u1>: (C =",8, ')\n', vec)
    print("------------------------------------------")
    x = dot(linalg.inv(M), b)
    statex = x/linalg.norm(x)
    print("|x>_true: \n", end="")
    print_state(statex)

main()