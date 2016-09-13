import numpy as np
from numpy import linalg

class Morses:

    def __init__(self, atomi, atomj, D, alpha, r0):
        self.atomi = atomi
        self.atomj = atomj
        self.D = D
        self.alpha = alpha
        self.r0 = r0
        self.r = None

class Bonds:

    def __init__(self, atomi, atomj, k, b0):
        self.atomi = atomi
        self.atomj = atomj
        self.k = k
        self.b0 = b0
        self.b = None

class Angles:

    def __init__(self, atomi, atomj, atomk, k, a0):
        self.atomi = atomi
        self.atomj = atomj
        self.atomk = atomk
        self.k = k
        self.a0 = a0
        self.a = None

class Dihedrals:

    def __init__(self, atomi, atomj, atomk, atoml, k, d0, n=None):
        self.atomi = atomi
        self.atomj = atomj
        self.atomk = atomk
        self.atoml = atoml
        self.k = k
        self.d0 = d0
        self.n = n
        self.d = None

def get_morse_potential_eta(atoms, morse):

    i = morse.atomi
    j = morse.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)

    if dij > morse.r0:
        exp = np.exp(-morse.alpha*(dij-morse.r0))
        eta = 1.0 - (1.0 - exp)**2
    else:
        eta = 1.0

    return eta

def get_morse_potential_value(atoms, morse):

    i = morse.atomi
    j = morse.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)

    exp = np.exp(-morse.alpha*(dij-morse.r0))

    v = morse.D*(1.0-exp)**2

    morse.r = dij

    return i, j, v

def get_morse_potential_gradient(atoms, morse):

    Mx=np.array([[1, 0, 0, -1, 0, 0],
                 [0, 1, 0, 0, -1, 0],
                 [0, 0, 1, 0, 0, -1]])

    i = morse.atomi
    j = morse.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij

    exp = np.exp(-morse.alpha*(dij-morse.r0))

    gr = 2.0*morse.D*morse.alpha*exp*(1.0-exp)*eij

    gx = np.dot(Mx.T, gr)

    morse.r = dij

    return i, j, gx

def get_morse_potential_hessian(atoms, morse, spectral=False):

    Mx=np.array([[1, 0, 0, -1, 0, 0],
                 [0, 1, 0, 0, -1, 0],
                 [0, 0, 1, 0, 0, -1]])

    i = morse.atomi
    j = morse.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij

    Pij = np.tensordot(eij,eij,axes=0)
    Qij = np.eye(3)-Pij

    exp = np.exp(-morse.alpha*(dij-morse.r0))

    Hr = 2.0*morse.D*morse.alpha*exp*(morse.alpha*(2.0*exp-1.0)*Pij + (1.0-exp)/dij*Qij)

    Hx = np.dot(Mx.T, np.dot(Hr, Mx))

    if spectral:
        eigvals, eigvecs = linalg.eigh(Hx)
        D = np.diag(np.abs(eigvals))
        U = eigvecs
        Hx = np.dot(U,np.dot(D,np.transpose(U)))

    morse.r = dij

    return i, j, Hx

def get_morse_potential_reduced_hessian(atoms, morse):

    Mx=np.array([[1, 0, 0, -1, 0, 0],
                 [0, 1, 0, 0, -1, 0],
                 [0, 0, 1, 0, 0, -1]])

    i = morse.atomi
    j = morse.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij

    Pij = np.tensordot(eij,eij,axes=0)

    exp = np.exp(-morse.alpha*(dij-morse.r0))

    Hr = np.abs(2.0*morse.D*morse.alpha**2*exp*(2.0*exp-1.0))*Pij

    Hx = np.dot(Mx.T, np.dot(Hr, Mx))

    morse.r = dij

    return i, j, Hx

def get_bond_potential_value(atoms, bond):

    i = bond.atomi
    j = bond.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)

    v = 0.5*bond.k*(dij-bond.b0)**2

    bond.b = dij

    return i, j, v

def get_bond_potential_gradient(atoms, bond):

    Bx=np.array([[1, 0, 0, -1, 0, 0],
                 [0, 1, 0, 0, -1, 0],
                 [0, 0, 1, 0, 0, -1]])

    i = bond.atomi
    j = bond.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij

    gr = bond.k*(dij-bond.b0)*eij

    gx = np.dot(Bx.T, gr)

    bond.b = dij

    return i, j, gx

def get_bond_potential_hessian(atoms, bond, morses=None, spectral=False):

    Bx=np.array([[1, 0, 0, -1, 0, 0],
                 [0, 1, 0, 0, -1, 0],
                 [0, 0, 1, 0, 0, -1]])

    i = bond.atomi
    j = bond.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij

    Pij = np.tensordot(eij,eij,axes=0)
    Qij = np.eye(3)-Pij

    Hr = bond.k*Pij+bond.k*(dij-bond.b0)/dij*Qij

    if morses is not None:
        for m in range(len(morses)):
            if morses[m].atomi == i or morses[m].atomi == j:
                Hr *= get_morse_potential_eta(atoms, morses[m])
            elif morses[m].atomj == i or morses[m].atomj == j:
                Hr *= get_morse_potential_eta(atoms, morses[m])

    Hx = np.dot(Bx.T, np.dot(Hr, Bx))

    if spectral:
        eigvals, eigvecs = linalg.eigh(Hx)
        D = np.diag(np.abs(eigvals))
        U = eigvecs
        Hx = np.dot(U,np.dot(D,np.transpose(U)))

    bond.b = dij

    return i, j, Hx

def get_bond_potential_reduced_hessian(atoms, bond, morses=None):

    Bx=np.array([[1, 0, 0, -1, 0, 0],
                 [0, 1, 0, 0, -1, 0],
                 [0, 0, 1, 0, 0, -1]])

    i = bond.atomi
    j = bond.atomj

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij

    Pij = np.tensordot(eij,eij,axes=0)

    Hr = bond.k*Pij

    if morses is not None:
        for m in range(len(morses)):
            if morses[m].atomi == i or morses[m].atomi == j:
                Hr *= get_morse_potential_eta(atoms, morses[m])
            elif morses[m].atomj == i or morses[m].atomj == j:
                Hr *= get_morse_potential_eta(atoms, morses[m])

    Hx = np.dot(Bx.T, np.dot(Hr, Bx))

    bond.b = dij

    return i, j, Hx

def get_bond_potential_reduced_hessian_test(atoms, bond):

    i, j, v = get_bond_potential_value(atoms, bond)
    i, j, gx = get_bond_potential_gradient(atoms, bond)

    Hx = np.tensordot(gx,gx,axes=0)/v/2.0

    return i, j, Hx

def get_angle_potential_value(atoms, angle):

    i = angle.atomi
    j = angle.atomj
    k = angle.atomk

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij
    rkj = rel_pos_pbc(atoms, k, j)
    dkj = linalg.norm(rkj)
    ekj = rkj/dkj
    eijekj = np.dot(eij, ekj)

    a = np.arccos(eijekj)
    da = a-angle.a0
    da = da - np.around(da / np.pi) * np.pi

    v = 0.5*angle.k*da**2

    angle.a = a

    return i, j, k, v

def get_angle_potential_gradient(atoms, angle):

    Ax=np.array([[1, 0, 0, -1, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, -1, 0, 0, 0, 0],
                 [0, 0, 1, 0, 0, -1, 0, 0, 0],
                 [0, 0, 0, -1, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, -1, 0, 0, 1, 0],
                 [0, 0, 0, 0, 0, -1, 0, 0, 1]])

    i = angle.atomi
    j = angle.atomj
    k = angle.atomk

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    eij = rij/dij
    rkj = rel_pos_pbc(atoms, k, j)
    dkj = linalg.norm(rkj)
    ekj = rkj/dkj
    eijekj = np.dot(eij, ekj)

    a = np.arccos(eijekj)
    da = a-angle.a0
    da = da - np.around(da / np.pi) * np.pi
    sina = np.sin(a)

    Pij = np.tensordot(eij,eij,axes=0)
    Qij = np.eye(3)-Pij
    Pkj = np.tensordot(ekj,ekj,axes=0)
    Qkj = np.eye(3)-Pkj

    gr = np.zeros(6)
    gr[0:3] = -angle.k*da/sina/dij*np.dot(Qij,ekj)
    gr[3:6] = -angle.k*da/sina/dkj*np.dot(Qkj,eij)

    gx = np.dot(Ax.T, gr)

    angle.a = a

    return i, j, k, gx

def get_angle_potential_hessian(atoms, angle, morses=None, spectral=False):

    Ax=np.array([[1, 0, 0, -1, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, -1, 0, 0, 0, 0],
                 [0, 0, 1, 0, 0, -1, 0, 0, 0],
                 [0, 0, 0, -1, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, -1, 0, 0, 1, 0],
                 [0, 0, 0, 0, 0, -1, 0, 0, 1]])

    i = angle.atomi
    j = angle.atomj
    k = angle.atomk

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    dij2 = dij*dij
    eij = rij/dij
    rkj = rel_pos_pbc(atoms, k, j)
    dkj = linalg.norm(rkj)
    dkj2 = dkj*dkj
    ekj = rkj/dkj
    dijdkj = dij*dkj
    eijekj = np.dot(eij, ekj)

    a = np.arccos(eijekj)
    da = a-angle.a0
    da = da - np.around(da / np.pi) * np.pi
    sina = np.sin(a)
    cosa = np.cos(a)
    ctga = cosa/sina

    Pij = np.tensordot(eij,eij,axes=0)
    Qij = np.eye(3)-Pij
    Pkj = np.tensordot(ekj,ekj,axes=0)
    Qkj = np.eye(3)-Pkj
    Pik = np.tensordot(eij,ekj,axes=0)
    Pki = np.tensordot(ekj,eij,axes=0)
    P = np.eye(3)*eijekj

    QijPkjQij = np.dot(Qij, np.dot(Pkj, Qij))
    QijPkiQkj = np.dot(Qij, np.dot(Pki, Qkj))
    QkjPijQkj = np.dot(Qkj, np.dot(Pij, Qkj))

    Hr = np.zeros((6,6))
    Hr[0:3,0:3] = angle.k*(QijPkjQij/sina + da*(-ctga*QijPkjQij/sina+np.dot(Qij, Pki)-np.dot(Pij, Pki)*2.0+(Pik+P)))/sina/dij2
    Hr[0:3,3:6] = angle.k*(QijPkiQkj/sina + da*(-ctga*QijPkiQkj/sina-np.dot(Qij, Qkj)))/sina/dijdkj
    Hr[3:6,0:3] = Hr[0:3,3:6].T
    Hr[3:6,3:6] = angle.k*(QkjPijQkj/sina + da*(-ctga*QkjPijQkj/sina+np.dot(Qkj, Pik)-np.dot(Pkj, Pik)*2.0+(Pki+P)))/sina/dkj2

    if morses is not None:
        for m in range(len(morses)):
            if morses[m].atomi == i or morses[m].atomi == j or morses[m].atomi == k:
                Hr *= get_morse_potential_eta(atoms, morses[m])
            elif morses[m].atomj == i or morses[m].atomj == j or morses[m].atomj == k:
                Hr *= get_morse_potential_eta(atoms, morses[m])

    Hx=np.dot(Ax.T, np.dot(Hr, Ax))

    if spectral:
        eigvals, eigvecs = linalg.eigh(Hx)
        D = np.diag(np.abs(eigvals))
        U = eigvecs
        Hx = np.dot(U,np.dot(D,np.transpose(U)))

    angle.a = a

    return i, j, k, Hx

def get_angle_potential_reduced_hessian(atoms, angle, morses=None):

    Ax=np.array([[1, 0, 0, -1, 0, 0, 0, 0, 0],
                 [0, 1, 0, 0, -1, 0, 0, 0, 0],
                 [0, 0, 1, 0, 0, -1, 0, 0, 0],
                 [0, 0, 0, -1, 0, 0, 1, 0, 0],
                 [0, 0, 0, 0, -1, 0, 0, 1, 0],
                 [0, 0, 0, 0, 0, -1, 0, 0, 1]])

    i = angle.atomi
    j = angle.atomj
    k = angle.atomk

    rij = rel_pos_pbc(atoms, i, j)
    dij = linalg.norm(rij)
    dij2 = dij*dij
    eij = rij/dij
    rkj = rel_pos_pbc(atoms, k, j)
    dkj = linalg.norm(rkj)
    dkj2 = dkj*dkj
    ekj = rkj/dkj
    dijdkj = dij*dkj
    eijekj = np.dot(eij, ekj)

    a = np.arccos(eijekj)
    sina = np.sin(a)
    sina2 = sina*sina

    Pij = np.tensordot(eij,eij,axes=0)
    Qij = np.eye(3)-Pij
    Pkj = np.tensordot(ekj,ekj,axes=0)
    Qkj = np.eye(3)-Pkj
    Pki = np.tensordot(ekj,eij,axes=0)

    Hr = np.zeros((6,6))
    Hr[0:3,0:3] = np.dot(Qij, np.dot(Pkj, Qij))/dij2
    Hr[0:3,3:6] = np.dot(Qij, np.dot(Pki, Qkj))/dijdkj
    Hr[3:6,0:3] = Hr[0:3,3:6].T
    Hr[3:6,3:6] = np.dot(Qkj, np.dot(Pij, Qkj))/dkj2

    Hr = Hr*angle.k/sina2

    if morses is not None:
        for m in range(len(morses)):
            if morses[m].atomi == i or morses[m].atomi == j or morses[m].atomi == k:
                Hr *= get_morse_potential_eta(atoms, morses[m])
            elif morses[m].atomj == i or morses[m].atomj == j or morses[m].atomj == k:
                Hr *= get_morse_potential_eta(atoms, morses[m])

    Hx=np.dot(Ax.T, np.dot(Hr, Ax))

    angle.a = a

    return i, j, k, Hx

def get_angle_potential_reduced_hessian_test(atoms, angle):

    i, j, k, v = get_angle_potential_value(atoms, angle)
    i, j, k, gx = get_angle_potential_gradient(atoms, angle)

    Hx = np.tensordot(gx,gx,axes=0)/v/2.0

    return i, j, k, Hx

def get_dihedral_potential_value(atoms, dihedral):

    i = dihedral.atomi
    j = dihedral.atomj
    k = dihedral.atomk
    l = dihedral.atoml

    rij = rel_pos_pbc(atoms, i, j)
    rkj = rel_pos_pbc(atoms, k, j)
    rkl = rel_pos_pbc(atoms, k, l)

    rmj = np.cross(rij, rkj)
    dmj = linalg.norm(rmj)
    emj = rmj/dmj
    rnk = np.cross(rkj, rkl)
    dnk = linalg.norm(rnk)
    enk = rnk/dnk
    emjenk = np.dot(emj, enk)

    d = np.sign(np.dot(rkj, np.cross(rmj, rnk)))*np.arccos(emjenk)
    dd = d-dihedral.d0
    dd = dd - np.around(dd / np.pi / 2.0) * np.pi * 2.0

    if dihedral.n is None:
        v = 0.5*dihedral.k*dd**2
    else:
        v = dihedral.k*(1.0 - np.cos(dihedral.n*d - dihedral.d0))

    dihedral.d = d

    return i, j, k, l, v

def get_dihedral_potential_gradient(atoms, dihedral):

    i = dihedral.atomi
    j = dihedral.atomj
    k = dihedral.atomk
    l = dihedral.atoml

    rij = rel_pos_pbc(atoms, i, j)
    rkj = rel_pos_pbc(atoms, k, j)
    dkj = linalg.norm(rkj)
    dkj2 = dkj*dkj
    rkl = rel_pos_pbc(atoms, k, l)

    rijrkj = np.dot(rij, rkj)
    rkjrkl = np.dot(rkj, rkl)

    rmj = np.cross(rij, rkj)
    dmj = linalg.norm(rmj)
    dmj2 = dmj*dmj
    emj = rmj/dmj
    rnk = np.cross(rkj, rkl)
    dnk = linalg.norm(rnk)
    dnk2 = dnk*dnk
    enk = rnk/dnk
    emjenk = np.dot(emj, enk)

    d = np.sign(np.dot(rkj, np.cross(rmj, rnk)))*np.arccos(emjenk)
    dd = d-dihedral.d0
    dd = dd - np.around(dd / np.pi / 2.0) * np.pi * 2.0

    dddri = dkj/dmj2*rmj
    dddrl = -dkj/dnk2*rnk

    gx = np.zeros(12)

    gx[0:3] = dddri
    gx[3:6] = (rijrkj/dkj2-1.0)*dddri-rkjrkl/dkj2*dddrl
    gx[6:9] = (rkjrkl/dkj2-1.0)*dddrl-rijrkj/dkj2*dddri
    gx[9:12] = dddrl

    if dihedral.n is None:
        gx *= dihedral.k*dd
    else:
        gx *= dihedral.k*dihedral.n*np.sin(dihedral.n*d - dihedral.d0)

    dihedral.d = d

    return i, j, k, l, gx

def get_dihedral_potential_hessian(atoms, dihedral, morses=None, spectral=False):

    eps = 0.000001

    i,j,k,l,g = get_dihedral_potential_gradient(atoms, dihedral)

    Hx = np.zeros((12,12))

    dihedral_eps = Dihedrals(dihedral.atomi, dihedral.atomj, dihedral.atomk, dihedral.atoml, dihedral.k, dihedral.d0, dihedral.n)
    indx = [3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2, 3*k, 3*k+1, 3*k+2, 3*l, 3*l+1, 3*l+2]
    for x in range(12):
        a = atoms.copy()
        positions = np.reshape(a.get_positions(),-1)
        positions[indx[x]] += eps
        a.set_positions(np.reshape(positions, (len(a),3)))
        i,j,k,l,geps = get_dihedral_potential_gradient(a, dihedral_eps)
        for y in range(12):
            Hx[x,y] += 0.5*(geps[y]-g[y])/eps
            Hx[y,x] += 0.5*(geps[y]-g[y])/eps

    if morses is not None:
        for m in range(len(morses)):
            if morses[m].atomi == i or morses[m].atomi == j or morses[m].atomi == k or morses[m].atomi == l:
                Hx *= get_morse_potential_eta(atoms, morses[m])
            elif morses[m].atomj == i or morses[m].atomj == j or morses[m].atomj == k or morses[m].atomj == l:
                Hx *= get_morse_potential_eta(atoms, morses[m])

    if spectral:
        eigvals, eigvecs = linalg.eigh(Hx)
        D = np.diag(np.abs(eigvals))
        U = eigvecs
        Hx = np.dot(U,np.dot(D,np.transpose(U)))

    return i, j, k, l, Hx

def get_dihedral_potential_reduced_hessian(atoms, dihedral, morses=None):

    i = dihedral.atomi
    j = dihedral.atomj
    k = dihedral.atomk
    l = dihedral.atoml

    rij = rel_pos_pbc(atoms, i, j)
    rkj = rel_pos_pbc(atoms, k, j)
    dkj = linalg.norm(rkj)
    dkj2 = dkj*dkj
    rkl = rel_pos_pbc(atoms, k, l)

    rijrkj = np.dot(rij, rkj)
    rkjrkl = np.dot(rkj, rkl)

    rmj = np.cross(rij, rkj)
    dmj = linalg.norm(rmj)
    dmj2 = dmj*dmj
    emj = rmj/dmj
    rnk = np.cross(rkj, rkl)
    dnk = linalg.norm(rnk)
    dnk2 = dnk*dnk
    enk = rnk/dnk
    emjenk = np.dot(emj, enk)

    d = np.sign(np.dot(rkj, np.cross(rmj, rnk)))*np.arccos(emjenk)

    dddri = dkj/dmj2*rmj
    dddrl = -dkj/dnk2*rnk

    gx = np.zeros(12)

    gx[0:3] = dddri
    gx[3:6] = (rijrkj/dkj2-1.0)*dddri-rkjrkl/dkj2*dddrl
    gx[6:9] = (rkjrkl/dkj2-1.0)*dddrl-rijrkj/dkj2*dddri
    gx[9:12] = dddrl

    if dihedral.n is None: 
        Hx = dihedral.k*np.tensordot(gx,gx,axes=0)
    else:
        Hx = np.abs(dihedral.k*dihedral.n**2*np.cos(dihedral.n*d - dihedral.d0))*np.tensordot(gx,gx,axes=0)

    if morses is not None:
        for m in range(len(morses)):
            if morses[m].atomi == i or morses[m].atomi == j or morses[m].atomi == k or morses[m].atomi == l:
                Hx *= get_morse_potential_eta(atoms, morses[m])
            elif morses[m].atomj == i or morses[m].atomj == j or morses[m].atomj == k or morses[m].atomj == l:
                Hx *= get_morse_potential_eta(atoms, morses[m])

    dihedral.d = d

    return i, j, k, l, Hx

def get_dihedral_potential_reduced_hessian_test(atoms, dihedral):

    i, j, k, l, gx = get_dihedral_potential_gradient(atoms, dihedral)

    if dihedral.n is None:
        i, j, k, l, v = get_dihedral_potential_value(atoms, dihedral)
        Hx = np.tensordot(gx,gx,axes=0)/v/2.0
    else:
        arg = dihedral.n*dihedral.d - dihedral.d0
        Hx = np.tensordot(gx,gx,axes=0)/dihedral.k/np.sin(arg)/np.sin(arg)*np.cos(arg)

    return i, j, k, l, Hx

def rel_pos_pbc(atoms, i, j):
    """
    Return difference between two atomic positions, correcting for jumps across PBC
    """
    d = atoms.get_positions()[i,:]-atoms.get_positions()[j,:]
    g = linalg.inv(atoms.get_cell().T)
    f = np.floor(np.dot(g, d.T) + 0.5)
    d -= np.dot(atoms.get_cell().T, f).T
    return d

def translational_vectors(atoms, mass_weighted=False):
    """
    Return normalised translational vectors
    """

    Tr = np.zeros((3*len(atoms),3))

    if mass_weighted:
        masses = atoms.get_masses()
    else:
        masses = np.ones(len(atoms))
    masses_sqrt = np.sqrt(masses)

    k=0
    for i in range(len(atoms)):
        for j in range(3):
           Tr[k,j] = masses_sqrt[i]
           k+=1

    for i in range(3):
        norm = np.sqrt(np.dot(Tr[:,i], Tr[:,i]))
        Tr[:,i] /= norm

    return Tr

def rotational_vectors(atoms, mass_weighted=False):
    """
    Return normalised rotational vectors
    """

    Rot = np.zeros((3*len(atoms),3))

    threshold = np.finfo(float).eps*10000.0

    if mass_weighted:
        masses = atoms.get_masses()
    else:
        masses = np.ones(len(atoms))
    masses_sqrt = np.sqrt(masses)

    com = np.zeros(3)
    for i in range(len(atoms)):
        com += masses[i] * atoms.get_positions()[i,:]
    com /= np.sum(masses)

    it = np.zeros((3,3))
    for i in range(len(atoms)):
        rpos = atoms.get_positions()[i,:] - com
        it[0,0] += masses[i] * (rpos[1]**2 + rpos[2]**2)
        it[1,1] += masses[i] * (rpos[0]**2 + rpos[2]**2)
        it[2,2] += masses[i] * (rpos[0]**2 + rpos[1]**2)
        it[0,1] -= masses[i] * (rpos[0]*rpos[1])
        it[0,2] -= masses[i] * (rpos[0]*rpos[2])
        it[1,2] -= masses[i] * (rpos[1]*rpos[2])
    it[1,0] = it[0,1]
    it[2,0] = it[0,2]
    it[2,1] = it[1,2]
    d, dit = linalg.eigh(it)

    for i in range(len(atoms)):
        rpos = atoms.get_positions()[i,:] - com
        cp = np.dot(np.transpose(dit), rpos)
        Rot[i*3,0] = masses_sqrt[i] * (cp[1]*dit[0,2]-cp[2]*dit[0,1])
        Rot[i*3+1,0] = masses_sqrt[i] * (cp[1]*dit[1,2]-cp[2]*dit[1,1])
        Rot[i*3+2,0] = masses_sqrt[i] * (cp[1]*dit[2,2]-cp[2]*dit[2,1])

        Rot[i*3,1] = masses_sqrt[i] * (cp[2]*dit[0,0]-cp[0]*dit[0,2])
        Rot[i*3+1,1] = masses_sqrt[i] * (cp[2]*dit[1,0]-cp[0]*dit[1,2])
        Rot[i*3+2,1] = masses_sqrt[i] * (cp[2]*dit[2,0]-cp[0]*dit[2,2])

        Rot[i*3,2] = masses_sqrt[i] * (cp[0]*dit[0,1]-cp[1]*dit[0,0])
        Rot[i*3+1,2] = masses_sqrt[i] * (cp[0]*dit[1,1]-cp[1]*dit[1,0])
        Rot[i*3+2,2] = masses_sqrt[i] * (cp[0]*dit[2,1]-cp[1]*dit[2,0])

    ndof = 3
    for i in range(3):
        norm = np.sqrt(np.dot(Rot[:,i], Rot[:,i]))
        if norm <= threshold:
            ndof -= 1
            continue
        Rot[:,i] /= norm
        if i < 2:
            for j in range(i+1):
                Rot[:,i+1] = Rot[:,i+1] - np.dot(Rot[:,i+1], Rot[:,j]) * Rot[:,j]

    return Rot[:,0:ndof]

def remove_tr_rot_vector(atoms, vecin, mass_weighted=False):

    Tr = translational_vectors(atoms, mass_weighted)
    Rot = rotational_vectors(atoms, mass_weighted)

    vecout = vecin

    for i in range(np.shape(Tr)[1]):
        norm = np.dot(vecout, Tr[:,i])
        vecout -= norm * Tr[:,i]

    for i in range(np.shape(Rot)[1]):
        norm = np.dot(vecout, Rot[:,i])
        vecout -= norm * Rot[:,i]

    return vecout
