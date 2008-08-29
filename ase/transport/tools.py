import numpy as npy 

def tri2full(M, UL='L'):
    """UP='L' => fill upper triangle from lower triangle
       such that M=M^d"""
    nbf = len(M)
    if UL == 'L':
        for i in range(nbf - 1):
            M[i, i:] = M[i:, i].conj()
    elif UL == 'U':
        for i in range(nbf - 1):
            M[i:, i] = M[i, i:].conj()

def dagger(matrix):
    return npy.conj(matrix.T)

def rotate_matrix(h, u):
    return npy.dot(u.T.conj(), npy.dot(h, u))

def get_subspace(matrix, index):
    """Get the subspace spanned by the basis function listed in index"""
    return npy.take(npy.take(matrix, index, axis=0), index, axis=1)   
    #return matrix[index, index]

def normalize_rot(c, s):
    """normalize column vectors so that <c[:,i]|s|c[:,i]> = 1 """
    for i in xrange(c.shape[0]):
        v = c[:, i]
        norm = 1.0 / npy.sqrt(npy.dot(v.conj(), npy.dot(s, v)))
        c[:, i] *= norm

def subdiagonalize(h_ii, s_ii, index_j):
    nb = h_ii.shape[0]
    nb_sub = len(index_j)
    h_sub_jj = get_subspace(h_ii, index_j)
    s_sub_jj = get_subspace(s_ii, index_j)
    e_j, v_jj = npy.linalg.eig(npy.linalg.solve(s_sub_jj, h_sub_jj))
    normalize_rot(v_jj, s_sub_jj) # normalize: <v_j|s|v_j> = 1
    permute_list = npy.argsort(e_j.real)
    e_j = npy.take(e_j, permute_list)
    v_jj = npy.take(v_jj, permute_list, axis=1)
    
    #setup transformation matrix
    c_ii = npy.identity(nb, complex)
    for i in xrange(nb_sub):
        for j in xrange(nb_sub):
            c_ii[index_j[i], index_j[j]] = v_jj[i, j]

    h1_ii = rotate_matrix(h_ii, c_ii)
    s1_ii = rotate_matrix(s_ii, c_ii)

    return h1_ii, s1_ii, c_ii, e_j

def cutcoupling(h, s, index_n):
    for i in index_n:
        s[:, i] = 0.0
        s[i, :] = 0.0
        s[i, i] = 1.0
        Ei = h[i, i]
        h[:, i] = 0.0
        h[i, :] = 0.0
        h[i, i] = Ei
