# reference: M.A.Nielsen and I.L.Chuang
#  "Quantum Computation and Quantum Information"
#   section 5.1

from math import pi

def QFT(qc, x, inv=False):
    """ Quantum Fourier Transform
    qc: QuantumCircuit object
    x: list(int)
      qubit indices
    inv: bool
      inverse transform or not
    output: qc
    """
    a = pi/2
    if inv: a = -a
    for i in range(len(x)-1,-1,-1):
        qc.h(x[i])
        for j in range(i):
            qc.cp(a/(1<<j), x[i-j-1], x[i])

    for i in range(len(x)>>1):
        qc.swap(x[i], x[-1-i])
