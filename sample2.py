# エンタングルメント生成回路　計算基底で測定

from qiskit import ClassicalRegister, QuantumRegister
from qiskit import QuantumCircuit, execute, BasicAer

qr = QuantumRegister(2)
cr = ClassicalRegister(2)
qc = QuantumCircuit(qr, cr)

qc.h(qr[0])
qc.cx(qr[0], qr[1])

qc.measure(qr[0], cr[0])
qc.measure(qr[1], cr[1])

backend = BasicAer.get_backend('qasm_simulator')
job = execute(qc, backend, shots=1024)
result = job.result()

print(result.get_counts(qc))
