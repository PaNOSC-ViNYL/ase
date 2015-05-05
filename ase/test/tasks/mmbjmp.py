"""Collection of Heusler systems.

From this paper:

  Markus Meinert

  Modified Becke-Johnson potential investigation of half-metallic Heusler compounds

  Phys. Rev. B 79, 085104 (2009) XXX

  http://arxiv.org/pdf/1210.7738v1.pdf
"""

# Data extracted with:
#     pdftotext -layout -f 2 -l 2 1210.7738v1.pdf - | "eaxp|mmBJ|Eg" \
#     | tr -d '\n' | sed -E 's/[ \t]+/ /g' \
#     | sed -n 's/NV/compound NV/p' \
#     | sed -n 's/Eg.*PBE/EgPBE/p' | sed -n 's/ mBJ.*Eg/ EgmBJ/p' \
#     | sed -e 's/^[ \t]*//' | sed -e 's/\s\+/,/g'
#
#     pdftotext -layout -f 2 -l 2 1210.7738v1.pdf - | grep -E "[0-9][0-9]*" \
#     | sed -n '/Co2 YZ Heusler compounds/,$p' \
#     | sed -n '/become half-metals in mBJLDA/q;p' \
#     | cut -c 62- | grep -v Heusler | sed '/^$/d' | sed -e 's/^[ \t]*//' \
#     | sed -n 's/2 /2/p' | sed 's/\xe2\x80\x94/-/g' | sed -e 's/\s\+/,/g'

from ase.atoms import string2symbols
from ase.lattice.spacegroup import crystal


class HeuslerMeinertCollection:
    data1 = """
    #compound,NV,aexp,mexp,mPBE,EgPBE,mmBJ,EgmBJ,c
    Co2TiAl,25,5.85,0.74,1.00,0.40,1.00,1.11,1.12
    Co2TiSn,26,6.08,1.96,2.00,0.47,2.00,1.16,1.17
    Co2VAl,26,5.72,1.95,2.00,0.36*,2.00,0.65,1.13
    Co2ZrSn,26,6.25,1.81,2.00,0.50,2.00,1.50,1.16
    Co2CrGa,27,5.81,3.01,3.04,0.39*,3.00,1.06,1.18
    Co2MnAl,28,5.75,4.04,4.03,0.61*,4.04,1.29*,1.14
    Co2MnSi,29,5.65,4.97,5.00,0.81,5.00,1.42,1.15
    Co2MnGe,29,5.75,4.93,5.00,0.57,5.00,1.49,1.19
    Co2MnSn,29,5.98,5.08,5.03,0.39*,5.04,1.36*,1.19
    Co2FeAl,29,5.73,4.96,4.99,0.06*,5.00,0.75,1.14
    Co2FeGa,29,5.74,5.04,5.02,0.02*,5.00,0.80,1.20
    Co2FeSi,30,5.64,5.97,5.47,0.11*,5.79,0.82*,1.16
    Co2FeGe,30,5.74,5.90,5.63,0.09*,5.98,0.90*,1.20
    Mn2VAl,22,5.92,1.94,2.00,0.28,2.00,0.48,1.09
    Mn2VGa,22,5.91,1.88,1.99,0.02*,2.00,0.27,1.15
    Fe2VAl,24,5.76,0.00,-,-,-,0.31,1.12
    Fe2VGa,24,5.78,0.00,-,-,-,0.39,1.17
    Fe2TiSn,24,6.09,0.00,-,-,-,0.69,1.16
    Ru2MnSb,28,6.20,4.40,4.03,0.28*,4.06,0.44*,1.19
    Ni2MnSn,31,6.05,4.05,4.03,-,4.17,-,1.19
    Cu2MnAl,32,5.95,3.60,3.51,-,3.50,-,1.13
    Cu2MnSn,33,6.17,4.11,3.86,-,3.91,-,1.19
    Cr2CoGa,24,5.80,0.35,0.08,0.19*,0.03,0.66*,1.17
    Mn2CoAl,26,5.84,1.95,2.00,0.43,2.00,0.68,1.12
    Mn2CoGe,27,5.80,2.99,3.00,0.36,3.00,0.76,1.18
    Fe2CoSi,29,5.65,4.99,4.96,-,5.00,0.57,1.16
    """

    data_ref = {}
    for l in data1.split():
        if 'compound' not in l:
            l1 = l.split(',')
            row = []
            for v in l1[1:]:
                if '*' in v:  # gaps which are above or below the Fermi energy
                    row.append(float(v[:-1]))
                elif v == '-':
                    row.append(float(0.0))
                else:
                    row.append(float(v))
            data_ref[l1[0]] = row

    data = data_ref.copy()
    names = [d.split(',')[0] for d in data1.split()][1:]
    labels = data1.split()[0].split(',')

    def __getitem__(self, name):
        d = self.data[name]
        # the index of label in labels less one
        # (compound is already as key in d)
        a = d[self.labels.index('aexp') - 1]
        if name in ['Cr2CoGa', 'Mn2CoAl', 'Mn2CoGe', 'Fe2CoSi']:
            # http://en.wikipedia.org/wiki/Space_group
            sg = '216'
        else:
            sg = '225'
        symbols = string2symbols(name)
        symbols.pop(0)
        b = crystal(symbols=symbols,
                    basis=[(1. / 4, 1. / 4, 1. / 4),
                           (0., 0., 0.),
                           (1. / 2, 1. / 2, 1. / 2)],
                    spacegroup=225,
                    cellpar=[a, a, a, 90, 90, 90],
                    primitive_cell=True)
        # set average moments on all atoms (add + 2.0)
        magmom = d[self.labels.index('mexp') - 1] + 2.0
        m = [magmom / len(b)] * len(b)
        # break spin symmetry between atoms no. 1 and 2
        m[1] = m[1] + m[2]
        m[2] = - m[2]
        b.set_initial_magnetic_moments(m)
        if 0:
            from ase.visualize import view
            view(b)

        return b

    def keys(self):
        return self.names

if __name__ == '__main__':
    HeuslerMeinertCollection()['Co2TiAl']
