from qiskit import ClassicalRegister, QuantumRegister
from qiskit import QuantumCircuit, execute, BasicAer
from numpy import *
import time
import sys

def hadamards(circ, qreg):
    for j in range(qreg.size):
        circ.h(qreg[j])

def ctrl_Upow(circ, ctrl_qreg, targ_qreg):
    for j in range(ctrl_qreg.size):
        for k in range(2**(ctrl_qreg.size-j-1)):
            ctrl_U(circ, ctrl_qreg[j], targ_qreg)

def ctrl_U(circ, ctrl_qbit, targ_qreg):    
    circ.cu1(2*pi*(11/16), ctrl_qbit, targ_qreg[0])
    circ.cx(ctrl_qbit, targ_qreg[0])
    circ.cu1(2*pi*(1/16), ctrl_qbit, targ_qreg[0])
    circ.cx(ctrl_qbit, targ_qreg[0])

def qft_INV(circ, qreg):
    for j in range(qreg.size):
        circ.h(qreg[j])
        for k in range(2, qreg.size-j+1):
            circ.cu1(-2*pi/(2**k), qreg[j+k-1], qreg[j])
    reverse(circ, qreg)

def reverse(circ, qreg):
    for j in range(qreg.size // 2):
        circ.swap(qreg[j], qreg[qreg.size-1-j])

def phase_estimation(circ, ctrl, targ):
    hadamards(circ, ctrl)
    ctrl_Upow(circ, ctrl, targ)
    qft_INV(circ, ctrl)

def main():
    m = int(sys.argv[1])
    qr = QuantumRegister(m, 'qr')
    cr = ClassicalRegister(m, 'cr')
    psi = QuantumRegister(1, 'psi')

    qc = QuantumCircuit(qr, psi, cr)

    psi_vec = array([1, 1]/sqrt(2) )
    qc.initialize(psi_vec, psi) 
    
    phase_estimation(qc, qr, psi)

    for j in range(m):
        qc.measure(qr[j], cr[j])

    backend = BasicAer.get_backend('qasm_simulator')
    job = execute(qc, backend, shots=1024)
    result = job.result()
    counts = result.get_counts(qc)
    print(counts)  

start_time = time.perf_counter()
main()
print(time.perf_counter() - start_time)