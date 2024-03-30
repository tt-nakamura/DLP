# reference:
#  N.D.Mermin "Quantum Computer Science"
#   section 3.7

from qiskit import QuantumCircuit,Aer
from math import gcd
from QFT import QFT
from intlib import ContFrac
from QComp import PowMulMod

def shor(n, shots=64):
    """ Shor's method to factor integer n
    n: int
      n>0, composite and square-free
    return d or None
      d: int, scalar
        factor of n, 1<d<n
      None is returned if factor is not found.
    """
    l = n.bit_length()
    m = l if n==15 else 2*l # Mermin, p82
    x = [i for i in range(m)]
    y = [i for i in range(m,m+l)]
    z = [i for i in range(m+l,m+2*l+2)] # ancilla
    sim = Aer.get_backend('qasm_simulator')
    a,M = 2, 1<<m
    while a<n:
        d = gcd(a,n)
        if d>1: return d
        qc = QuantumCircuit(z[-1]+1, m)
        qc.h(x)
        qc.x(y[0])
        PowMulMod(qc, y, a, x, n, z)
        QFT(qc, x)
        qc.measure(x,x)
    
        job = sim.run(qc, shots=shots)
        c = job.result().get_counts(qc) # little endian
        r = {ContFrac(int(c,2),M,n)[1] for c in c}

        for r in sorted(r):
            if pow(a,r,n) > 1: continue
            while r&1 == 0:
                r >>= 1
                d = gcd(n, pow(a,r,n)+1)
                if d>1 and d<n: return d

        a += 1

if __name__ == '__main__':
    print(shor(15))
