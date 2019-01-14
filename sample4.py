# エンタングルメント生成回路　量子回路の表示

from qiskit import ClassicalRegister, QuantumRegister
from qiskit import QuantumCircuit, execute, BasicAer

qr = QuantumRegister(2, 'my qubit')
qc = QuantumCircuit(qr)

qc.h(qr[0])
qc.cx(qr[0], qr[1])

print(qc)