import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from QComp import add,sub

a,l = 7,4
x = [i for i in range(l)]

qc = QuantumCircuit(l)
add(qc, x, a)

plt.subplot(211)
qc.draw('mpl', ax=plt.gca())

qc = QuantumCircuit(l)
sub(qc, x, a)

plt.subplot(212)
qc.draw('mpl', ax=plt.gca())

plt.tight_layout()
plt.savefig('fig1.png')