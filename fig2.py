import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit,Aer
from QFT import QFT
from QComp import PowMulMod

a,m = 3,7

l = m.bit_length()
qc = QuantumCircuit(4*l+2, l)
x = [i for i in range(2*l)]
y = [i for i in range(2*l,3*l)]
z = [i for i in range(3*l,4*l+2)] # ancilla
qc.h(x)
qc.x(y[0])
PowMulMod(qc, y, a, x, m, z)
QFT(qc, x)
qc.measure(y, [i for i in range(l)])
   
sim = Aer.get_backend('statevector_simulator')
res = sim.run(qc).result()
c = res.get_counts(qc)
s = res.get_statevector(qc)
i = int(list(c.keys())[0], 2)<<len(x)
j = i + (1<<len(x))

plt.figure(figsize=(5,3.75))
plt.plot(np.abs(s[i:j])**2, '.:')
plt.xlabel('$j$ in equation(11)')
plt.ylabel('probablity distribution')
plt.tight_layout()
plt.savefig('fig2.eps')
