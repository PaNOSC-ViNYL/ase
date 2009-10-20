from ase.atoms import Atoms
import numpy as np

def gcd(m,n):
    while n:
        m,n=n,m%n
    return m

def nanotube(n, m, length=1, bond=1.42, symbol='C', verbose=False):
    if n < m:
        m, n = n, m
    nk = 6000
    sq3 = np.sqrt(3.0)
    a = sq3 * bond
    l2 = n * n + m * m + n * m
    l = np.sqrt(l2)
    dt = a * l / np.pi

    nd = gcd(n ,m)
    if (n - m) % (3 * nd ) == 0:
        ndr = 3 * nd
    else:
        ndr = nd

    nr = (2 * m + n) / ndr
    ns = -(2 * n + m) / ndr
    nt2 = 3 * l2 / ndr / ndr
    nt = np.floor(np.sqrt(nt2))
    nn = 2 * l2 / ndr

    ichk = 0
    if nr == 0:
        n60 = 1
    else:
        n60 = nr * 4
    
    absn = abs(n60)
    nnp = []
    nnq = []
    for i in range(-absn, absn + 1):
        for j in range(-absn, absn + 1):
            j2 = nr * j - ns * i
            if j2 == 1:
                j1 = m * i - n * j
                if j1 > 0 and j1 < nn:
                    ichk += 1
                    nnp.append(i)
                    nnq.append(j)

    if ichk == 0:
        raise RuntimeError('not found p, q strange!!')
    if ichk >= 2:
        raise RuntimeError('more than 1 pair p, q strange!!')

    nnnp = nnp[0]
    nnnq = nnq[0]

    if verbose:   
        print 'the symmetry vector is', nnnp, nnnq

    lp = nnnp * nnnp + nnnq * nnnq + nnnp * nnnq
    r = a * np.sqrt(lp)
    c = a * l
    t = sq3 * c / ndr
    
    if 2 * nn > nk:
        raise RuntimeError('parameter nk is too small!')

    rs = c / (2.0 * np.pi)
    
    if verbose:
        print 'radius=', rs, t

    q1 = np.arctan((sq3 * m) / (2 * n + m))
    q2 = np.arctan((sq3 * nnnq) / (2 * nnnp + nnnq))
    q3 = q1 - q2

    q4 = 2.0 * np.pi / nn
    q5 = bond * np.cos((np.pi / 6.0) - q1) / c * 2.0 * np.pi

    h1 = abs(t) / abs(np.sin(q3))
    h2 = bond * np.sin((np.pi / 6.0) - q1)

    ii = 0
    x, y, z = [], [], []    
    for i in range(nn):
        x1, y1, z1 = 0, 0, 0
    
        k = np.floor(i * abs(r) / h1)
        x1 = rs * np.cos(i * q4)
        y1 = rs * np.sin(i * q4)
        z1 = (i * abs(r) - k * h1) * np.sin(q3)
        kk2 = abs(np.floor((z1 + 0.0001) / t))
        if z1 >= t - 0.0001:
            z1 -= t * kk2
        elif z1 < 0:
            z1 += t * kk2
        ii += 1

        x.append(x1)
        y.append(y1)
        z.append(z1)
        z3 = (i * abs(r) - k * h1) * np.sin(q3) - h2
        ii += 1

        if z3 >= 0 and z3 < t:
            x2 = rs * np.cos(i * q4 + q5)
            y2 = rs * np.sin(i * q4 + q5)
            z2 = (i * abs(r) - k * h1) * np.sin(q3) - h2
            x.append(x2)
            y.append(y2)
            z.append(z2)
        else:
            x2 = rs * np.cos(i * q4 + q5)
            y2 = rs * np.sin(i * q4 + q5)
            z2 = (i * abs(r) - (k + 1) * h1) * np.sin(q3) - h2
            kk = abs(np.floor(z2 / t))
            if z2 >= t - 0.0001:
                z2 -= t * kk
            elif z2 < 0:
                z2 += t * kk
            x.append(x2)
            y.append(y2)
            z.append(z2)  
        
    ntotal = 2 * nn
    X = []
    for i in range(ntotal):
        X.append([x[i], y[i], z[i]])

    if length > 1:
        xx = X[:]
        for mnp in range(2, length + 1):
            for i in range(len(xx)):
                X.append(xx[i][:2] + [xx[i][2] + (mnp - 1) * t])
                
    TransVec = t
    NumAtom = ntotal * length
    Diameter = rs * 2
    ChiralAngle = np.arctan((sq3 * n) / (2 * m + n)) / (np.pi * 180)
    
    cell = [Diameter * 2, Diameter * 2, length * t]
    atoms = Atoms(symbol + str(NumAtom), positions=X)
    atoms.center()
    if verbose:
        print 'translation vector =', TransVec
        print 'diameter = ', Diameter
        print 'chiral angle = ', ChiralAngle
    return atoms

def graphene_nanoribbon(n, m, type='zigzag', saturated=False, C_H=1.09,
              C_C=1.42, vacc=5., magnetic=None, initial_mag=1.12):
    #This function creates the coordinates for a graphene nanoribbon,
    #n is width, m is length
    
    b = np.sqrt(3) * C_C / 4
    arm_unit = Atoms('C4', pbc=(1,0,1), cell = [4 * b,  vacc,  3 * C_C])
    arm_unit.positions = [[0, 0, 0],
                          [b * 2, 0, C_C / 2.],
                          [b * 2, 0, 3 * C_C / 2.],
                          [0, 0, 2 * C_C]]
    zz_unit = Atoms('C2', pbc=(1,0,1), cell = [3 * C_C /2., vacc, b * 4])
    zz_unit.positions = [[0, 0, 0],
                         [C_C / 2., 0, b * 2]]
    atoms = Atoms()    
    tol = 1e-4    
    if type == 'zigzag':
        edge_index0 = np.arange(m) * 2 + 1
        edge_index1 = (n - 1) * m * 2 + np.arange(m) * 2 
        if magnetic:
            mms = np.zeros(m * n * 2)
            for i in edge_index0:
                mms[i] = initial_mag * factor
            for i in edge_index1:
                mms[i] = -initial_mag * factor
                
        for i in range(n):
            layer = zz_unit.repeat((1, 1, m))
            layer.positions[:, 0] -= 3 * C_C / 2 * i
            if i % 2 == 1:
                layer.positions[:, 2] += 2 * b
                layer[-1].position[2] -= b * 4 * m
            atoms += layer
        if magnetic:
            atoms.set_initial_magnetic_moments(mms)
        if saturated:
            H_atoms0 = Atoms('H' + str(m))
            H_atoms0.positions = atoms[edge_index0].positions
            H_atoms0.positions[:, 0] += C_H
            H_atoms1 = Atoms('H' + str(m))
            H_atoms1.positions = atoms[edge_index1].positions
            H_atoms1.positions[:, 0] -= C_H
            atoms += H_atoms0 + H_atoms1
        atoms.cell = [n * 3 * C_C / 2 + vacc, vacc, m * 4 * b]
    
    elif type == 'armchair':
        for i in range(n):
            layer = arm_unit.repeat((1, 1, m))
            layer.positions[:, 0] -= 4 * b * i
            atoms += layer
        atoms.cell = [b * 4 * n + vacc, vacc, 3 * C_C * m]
    atoms.center()
    return atoms

