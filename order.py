# refrence: M.A.Nielsan and I.L.Chuang
#  "Quantum Computation and Quantum Information"
#   section 5.3

from qiskit import QuantumCircuit,Aer
from QFT import QFT
from intlib import ContFrac
from QComp import PowMulMod

def order(a, m, shots=64):
    """ order of a mod m
    a,m: int, scalar
    return r or None
      r: int, scalar
        a^r == 1 (mod m) and r is least
      None is returned if order is not found.
    """
    l = m.bit_length()
    qc = QuantumCircuit(4*l+2, 2*l)
    x = [i for i in range(2*l)]
    y = [i for i in range(2*l,3*l)]
    z = [i for i in range(3*l,4*l+2)] # ancilla
    qc.h(x)
    qc.x(y[0])
    PowMulMod(qc, y, a, x, m, z)
    QFT(qc, x)
    qc.measure(x,x)
    
    sim = Aer.get_backend('qasm_simulator')
    job = sim.run(qc, shots=shots)
    c = job.result().get_counts(qc) # little endian
    d = 1<<len(x)
    r = {ContFrac(int(c,2),d,m)[1] for c in c}

    for r in sorted(r):
        if pow(a,r,m) == 1: return r

if __name__ == '__main__':
    for i in range(2,7):
        print(order(i,7))
