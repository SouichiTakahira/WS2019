# エンタングルメント生成回路　ユニタリ行列で取得

from qiskit import ClassicalRegister, QuantumRegister
from qiskit import QuantumCircuit, execute, BasicAer

qr = QuantumRegister(2)
qc = QuantumCircuit(qr)

qc.h(qr[0])
qc.cx(qr[0], qr[1])

backend = BasicAer.get_backend('unitary_simulator')
job = execute(qc, backend)
result = job.result()

print(result.get_unitary(qc))