import numpy as np

def get_band_gap(calc):
    """Calculates the band gap - direct and indirect along with the 
    relevant k-points. 
    Rerurns a list of [e_direct, e_indirect, (is ik), (ivs, ivk), (ics, ick)],
    where (is, ik) are the irreducible k-point/spin indeces of the direct gap 
    and (ivs, ivk) and (ics, ick) are the irreducible spin and k-point indices 
    for the indirect gap in the valence and conduction bands respectively"""
    
    kpts_kc = calc.get_ibz_k_points()
    Nk = len(kpts_kc)
    Ns = calc.get_number_of_spins()
    e_skn = np.array([[calc.get_eigenvalues(kpt=k, spin=s)
                       for k in range(Nk)]
                      for s in range(Ns)])
    assert len(e_skn[0,0]) > 1

    e_skn -= calc.get_fermi_level()
    ev_sk = []
    ec_sk = []
    for s in range(Ns):
        for k in range(Nk):
            ev_n = e_skn[s, k][e_skn[s, k] < 0]
            if len(ev_n) == 0:
                ev_sk.append(-np.inf)
            else:
                ev_sk.append(np.max(ev_n))
            ec_n = e_skn[s, k][e_skn[s, k] > 0]
            if len(ec_n) == 0:
                ec_sk.append(np.inf)
            else:
                ec_sk.append(np.min(ec_n))
    ev_sk, ec_sk = np.array(ev_sk), np.array(ec_sk)

    edir = np.min(ec_sk - ev_sk)
    ein = np.min(ec_sk) - np.max(ev_sk)

    kdir = np.argmin(ec_sk - ev_sk)
    sdir = kdir % Ns
    kdir = kdir % Nk

    kvin = np.argmax(ev_sk)
    svin = kvin % Ns
    kvin = kvin % Nk

    kcin = np.argmin(ec_sk)
    scin = kcin % Ns
    kcin = kcin % Nk

    print 
    print 'Direct gap: %.3f eV' % edir
    k_c = kpts_kc[kdir]
    k_c = (k_c[0], k_c[1], k_c[2])
    print 'Transition at spin   :   %d' % sdir
    print 'Transition at k-point:  [%.3f %.3f %.3f]' % k_c
    print 
    print 'Indirect gap: %.3f eV' % ein
    kv_c = kpts_kc[kvin]
    kc_c = kpts_kc[kcin]
    k_2c = (kv_c[0], kv_c[1], kv_c[2], kc_c[0], kc_c[1], kc_c[2])
    print 'Transition (v -> c) at spin   :   %d  ->  %d' % (svin, scin)
    print 'Transition (v -> c) at k-point: ',
    print '[%.3f %.3f %.3f]  ->  [%.3f %.3f %.3f]' % k_2c
    print 
    return [edir, ein, (sdir, kdir), (svin, kvin), (scin, kcin)]
