import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit,Aer
from QFT import QFT
from QComp import PowMulMod

a,b,m = 3,5,7

l = m.bit_length()
qc = QuantumCircuit(6*l+2, l)
x = [i for i in range(2*l)]
y = [i for i in range(2*l,4*l)]
z = [i for i in range(4*l,5*l)]
w = [i for i in range(5*l,6*l+2)] # ancilla
qc.h(x+y)
qc.x(z[0])
PowMulMod(qc, z, a, x, m, w)
PowMulMod(qc, z, b, y, m, w)
QFT(qc, x)
QFT(qc, y)
qc.measure(z, [i for i in range(l)])
   
sim = Aer.get_backend('statevector_simulator')
res = sim.run(qc).result()
c = res.get_counts(qc)
s = res.get_statevector(qc)
i = int(list(c.keys())[0], 2)<<(4*l)
j = i + (1<<(4*l))
s = s[i:j].reshape(1<<(2*l), 1<<(2*l))

plt.figure(figsize=(5,4.1))
plt.contourf(np.abs(s), cmap='Reds')
plt.colorbar()
plt.xlabel('$j$ in equation(17)')
plt.ylabel('$k$ in equation(17)')
plt.axis('equal')
plt.tight_layout()
plt.savefig('fig3.eps')
