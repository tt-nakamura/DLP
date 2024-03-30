# refrence: M.A.Nielsan and I.L.Chuang
#  "Quantum Computation and Quantum Information"
#   section 5.4

from qiskit import QuantumCircuit,Aer
from QFT import QFT
from intlib import ContFrac,InvMod
from QComp import PowMulMod

def DLP(a, b, m, shots=64):
    """ discrete log problem to solve a^x == b (mod m)
    a,b,m: int, scalar
    shots: int, scalar
      number of measurement times
    return x or None
      x: int, scalar
        solution of DLP
      None is returned if solution is not found
    """
    l = m.bit_length()
    qc = QuantumCircuit(6*l+2, 4*l)
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
    qc.measure(x,x)
    qc.measure(y,y)
    
    sim = Aer.get_backend('qasm_simulator')
    job = sim.run(qc, shots=shots)
    c = job.result().get_counts(qc)
    q,d = set(), 1<<len(x)
    for c in c: # little endian
        j = int(c[len(x):], base=2) # a
        k = int(c[:len(x)], base=2) # b
        u,v = ContFrac(j,d,m)
        s,t = ContFrac(k,d,m)
        s,t = divmod(s*v, t)
        if u==0 or t: continue
        q.add(s*InvMod(u,v) % v)

    for q in sorted(q):
        if pow(a,q,m) == b: return q

if __name__ == '__main__':
    for i in range(2,7):
        print(DLP(3,i,7))
