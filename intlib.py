def ContFrac(a,b,m):
    """ continued fraction expansion of a/b
    a,b,m: int, scalar
      assume a>0, b>0, m>0
    return u,v: int, scalar
      u/v = rational approximation of a/b
        with largest v<m, |u/v - a/b| < 1/v^2,
        u and v are coprime, and v>0.
    """
    s,t,u,v = 0,1,1,0
    while b and v<m:
        q,r = divmod(a,b)
        a,b,s,t,u,v = b,r,u,v,s+u*q,t+v*q

    if v<m: return u,v
    else: return s,t

def InvMod(a,m):
    """ a^{-1} mod m
    a,m: int, scalar
      assume a and m are coprime
    return b: int, scalar
      inverse of a mod m, i.e., ab==1 (mod m)
    """
    s,u = 1,0
    while m:
        q,r = divmod(a,m)
        a,m,s,u = m,r,u,s-q*u
    if a!=1:
        raise RuntimeError("a,m not coprime")
    return s
