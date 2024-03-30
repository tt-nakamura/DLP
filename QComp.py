# reference:
#  V.Vedral, A.Barenco and A.Ekert
#  "Quantum Networks for Elementary Arithmetic Operations"
#   Physical Review A 54 (1996) 147

from intlib import InvMod

def add(qc, x, a, c=[]):
    """ x += a
    qc: QuantumCircuit object
    x: list(int)
      qubit indices
    a: int, scalar
      assume len(x) > a.bit_length()
    c: list(int)
      control qubit indices
    output: qc
    """
    for i in range(a.bit_length()):
        if a&(1<<i):
            for j in range(len(x)-1, i-1, -1):
                k = x[i:j] + c
                if k: qc.mcx(k, x[j])
                else: qc.x(x[j])

def sub(qc, x, a, c=[]):
    """ x -= a (same as add) """
    for i in range(a.bit_length()-1,-1,-1):
        if a&(1<<i):
            for j in range(i, len(x)):
                k = x[i:j] + c
                if k: qc.mcx(k, x[j])
                else: qc.x(x[j]) 

def AddMod(qc, x, a, m, s, c=[]):
    """ x = x + a (mod m)
    qc: QuantumCircuit object
    x: list(int)
      qubit indices
    a,m: int, scalar
      assume len(x) > m.bit_length() and 0<a<m
    s: int, scalar
      index of ancilla qubit
    c: list(int)
      control qubit indices
    output: qc
    ref: Vedral, Barenco and Ekert, fig4
    """
    k = x[-1:] + c
    add(qc, x, a, c)
    sub(qc, x, m, c)
    qc.mcx(k, s)
    add(qc, x, m, [s] + c)
    sub(qc, x, a, c)
    qc.mcx(k, s)
    if c: qc.mcx(c, s)
    else: qc.x(s)
    add(qc, x, a, c)
   

def SubMod(qc, x, a, m, s, c=[]):
    """ x = x + a (mod m) (same as AddMod) """
    k = x[-1:] + c
    sub(qc, x, a, c)
    if c: qc.mcx(c, s)
    else: qc.x(s)
    qc.mcx(k, s)
    add(qc, x, a, c)
    sub(qc, x, m, [s] + c)
    qc.mcx(k, s)
    add(qc, x, m, c)
    sub(qc, x, a, c)

def MulAddMod(qc, y, x, a, m, s, c=[]):
    """ y = y + x*a (mod m) (x unchanged)
    qc: QuantumCircuit object
    x,y: list(int)
      qubit indices
    a,m: int, scalar
      assume len(x) == m.bit_length() and
             len(y) == len(x) + 1 and 0<a<m
    s: int, scalar
      index of ancilla qubit
    c: list(int)
      control qubit indices    
    output: qc
    ref: Vedral, Barenco and Ekert, fig5
    """
    for i in x:
        AddMod(qc, y, a, m, s, [i] + c)
        a = (a<<1) % m

def MulSubMod(qc, y, x, a, m, s, c=[]):
    """ y = y - x*a (mod m) (same as MulAddMod) """
    for i in x:
        SubMod(qc, y, a, m, s, [i] + c)
        a = (a<<1) % m

def PowMulMod(qc, y, a, x, m, s):
    """ y = y * a**x (mod m) (x unchanged)
    qc: QuantumCircuit object
    x,y: list(int)
      qubit indices
    a,m: int, scalar
      assume len(x) > m.bit_length() and
             len(y) == m.bit_length() and 0<a<m
    s: list(int)
      indices of ancilla qubit
      assume len(s) == len(y) + 2
    output: qc
    ref: Vedral, Barenco and Ekert, fig6
    """
    s,z = s[0],s[1:]
    for i in x:
        MulAddMod(qc, z, y, a, m, s, [i])
        for j,k in zip(y,z): qc.cswap(i,j,k)
        b = InvMod(a,m)
        MulSubMod(qc, z, y, b, m, s, [i])
        a = (a*a) % m
