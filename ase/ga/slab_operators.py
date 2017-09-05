"""Operators that work on slabs.
Allowed compositions are respected."""
import random
import numpy as np

from ase.ga.offspring_creator import OffspringCreator
from ase.ga.element_mutations import get_row_column


def permute2(atoms):
    i1 = random.choice(range(len(atoms)))
    sym1 = atoms[i1].symbol
    i2 = random.choice([a.index for a in atoms if a.symbol != sym1])
    atoms[i1].symbol = atoms[i2].symbol
    atoms[i2].symbol = sym1


def get_add_remove_lists(**kwargs):
    to_add, to_rem = [], []
    for s, amount in kwargs.items():
        if amount > 0:
            to_add.extend([s] * amount)
        elif amount < 0:
            to_rem.extend([s] * abs(amount))
    return to_add, to_rem


def same_layer_comp(atoms):
    unique_syms, comp = np.unique(sorted(atoms.get_chemical_symbols()),
                                  return_counts=True)
    l = get_layer_comps(atoms)
    sym_dict = dict((s, np.array(c) / len(l))
                    for s, c in zip(unique_syms, comp))
    for la in l:
        correct_by = sym_dict.copy()
        lcomp = dict(
            zip(*np.unique([atoms[i].symbol for i in la], return_counts=True)))
        for s, num in lcomp.items():
            correct_by[s] -= num
        to_add, to_rem = get_add_remove_lists(**correct_by)
        for add, rem in zip(to_add, to_rem):
            ai = random.choice([i for i in la if atoms[i].symbol == rem])
            atoms[ai].symbol = add


def get_layer_comps(atoms, eps=1e-2):
    lc = []
    old_z = np.inf
    for z, ind in sorted([(a.z, a.index) for a in atoms]):
        if abs(old_z - z) < eps:
            lc[-1].append(ind)
        else:
            lc.append([ind])
        old_z = z

    return lc


class SlabOperator(OffspringCreator):
    def __init__(self, verbose=False, num_muts=1,
                 allowed_compositions=None,
                 same_layer_composition=True, element_pool=None):
        OffspringCreator.__init__(self, verbose, num_muts=num_muts)

        self.same_layer_composition = same_layer_composition
        self.allowed_compositions = allowed_compositions
        self.element_pool = element_pool

    def get_symbols_to_use(self, syms):
        """Get the symbols to use for the offspring candidate. The returned
        list of symbols will respect self.allowed_compositions"""
        if self.allowed_compositions is None:
            return syms

        unique_syms, counts = np.unique(syms, return_counts=True)
        comp, unique_syms = zip(*sorted(zip(counts, unique_syms),
                                        reverse=True))

        for cc in self.allowed_compositions:
            comp += (0,) * (len(cc) - len(comp))
            if comp == tuple(sorted(cc)):
                return syms

        comp_diff = self.get_closest_composition_diff(comp)
        to_add, to_rem = get_add_remove_lists(
            **dict(zip(unique_syms, comp_diff)))
        for add, rem in zip(to_add, to_rem):
            tbc = [i for i in range(len(syms)) if syms[i] == rem]
            ai = random.choice(tbc)
            syms[ai] = add
        return syms

    def get_closest_composition_diff(self, c):
        comp = np.array(c)
        mindiff = 1e10
        allowed_list = list(self.allowed_compositions)
        random.shuffle(allowed_list)
        for ac in allowed_list:
            difflen = len(comp) - len(ac)
            if difflen > 0:
                ac += (0,) * difflen
            diff = np.array(ac) - comp
            numdiff = sum([abs(i) for i in diff])
            if numdiff < mindiff:
                mindiff = numdiff
                ccdiff = diff
        return ccdiff

    def get_possible_mutations(self, a):
        unique_syms, comp = np.unique(sorted(a.get_chemical_symbols()),
                                      return_counts=True)
        min_num = min([i for i in np.ravel(list(self.allowed_compositions))
                       if i > 0])
        muts = set()
        for i, n in enumerate(comp):
            if n != 0:
                muts.add((unique_syms[i], n))
            if n % min_num >= 0:
                for j in range(1, n // min_num):
                    muts.add((unique_syms[i], min_num * j))
        return list(muts)

    def finalize_individual(self, indi):
        atoms_string = ''.join(indi.get_chemical_symbols())
        indi.info['key_value_pairs']['atoms_string'] = atoms_string
        return OffspringCreator.finalize_individual(self, indi)


class CutSpliceSlabCrossover(SlabOperator):
    def __init__(self, verbose=False, num_muts=1, tries=1000, min_ratio=0.25,
                 allowed_compositions=None, same_layer_composition=True):
        SlabOperator.__init__(self, verbose, num_muts,
                              allowed_compositions, same_layer_composition)

        self.tries = tries
        self.min_ratio = min_ratio
        self.descriptor = 'CutSpliceSlabCrossover'

    def get_new_individual(self, parents):
        f, m = parents

        indi = self.initialize_individual(f, f)
        indi.info['data']['parents'] = [i.info['confid'] for i in parents]

        fsp = f.get_scaled_positions()
        ma = np.max(fsp.transpose(), axis=1)
        mi = np.min(fsp.transpose(), axis=1)

        for _ in range(self.tries):
            # Find center point of cut
            rv = [random.random() for _ in range(3)]  # random vector
            midpoint = (ma - mi) * rv + mi

            # Determine cut plane
            theta = random.random() * 2 * np.pi  # 0,2pi
            phi = random.random() * np.pi  # 0,pi
            e = np.array((np.sin(phi) * np.cos(theta),
                          np.sin(theta) * np.sin(phi),
                          np.cos(phi)))

            # Cut structures
            d2fsp = np.dot(fsp - midpoint, e)
            fpart = d2fsp > 0
            ratio = float(np.count_nonzero(fpart)) / len(f)
            if ratio < self.min_ratio or ratio > 1 - self.min_ratio:
                continue
            syms = np.where(fpart, f.get_chemical_symbols(),
                            m.get_chemical_symbols())
            syms = self.get_symbols_to_use(syms)
            indi.set_chemical_symbols(syms)
            if self.same_layer_composition:
                same_layer_comp(indi)
            break

        # # The cell is set as the average of the two parents
        # # This may not be correct if the cut did not satisfy
        # # the compositional conditions
        # indi.set_cell((ratio * f.cell + (1 - ratio) * m.cell),
        #               scale_atoms=True)

        parent_message = ': Parents {0} {1}'.format(f.info['confid'],
                                                    m.info['confid'])
        return (self.finalize_individual(indi),
                self.descriptor + parent_message)

# Mutations: Random, MoveUp/Down/Left/Right, six or all elements


class RandomSlabMutation(SlabOperator):
    def __init__(self, verbose=False, num_muts=1, element_pool=None,
                 allowed_compositions=None, same_layer_composition=True):
        SlabOperator.__init__(self, verbose, num_muts,
                              allowed_compositions,
                              element_pool=element_pool)

        self.descriptor = 'RandomSlabMutation'

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f, f)
        indi.info['data']['parents'] = [i.info['confid'] for i in parents]

        # Do the operation
        muts = self.get_possible_mutations(f)
        tbmf = random.choice(muts)  # to be mutated from
        # tbmt = to be mutated to
        tbmt = random.choice(np.setdiff1d(self.element_pool,
                                          np.array([tbmf[0]])))
        # itbm = index to be mutated
        itbm = random.sample([a.index for a in indi if a.symbol == tbmf[0]],
                             tbmf[1])
        for i in itbm:
            indi[i].symbol = tbmt
        syms = self.get_symbols_to_use(indi.get_chemical_symbols())
        indi.set_chemical_symbols(syms)
        if self.same_layer_composition:
            same_layer_comp(indi)

        parent_message = ': Parent {0}'.format(f.info['confid'])
        return (self.finalize_individual(indi),
                self.descriptor + parent_message)


class NeighborhoodSlabMutation(SlabOperator):
    def __init__(self, verbose=False, num_muts=1, element_pool=None,
                 allowed_compositions=None, same_layer_composition=True):
        SlabOperator.__init__(self, verbose, num_muts,
                              allowed_compositions,
                              element_pool=element_pool)

        self.descriptor = 'NeighborhoodSlabMutation'

    def get_new_individual(self, parents):
        f = parents[0]

        indi = self.initialize_individual(f, f)
        indi.info['data']['parents'] = [i.info['confid'] for i in parents]

        # Do the operation
        muts = self.get_possible_mutations(f)
        tbmf = random.choice(muts)  # to be mutated from
        poss_muts = np.setdiff1d(self.element_pool, np.array([tbmf[0]]))
        rc = get_row_column(tbmf[0])
        random.shuffle(poss_muts)
        min_dist = 1e44
        for m in poss_muts:
            dist = sum(np.abs(np.array(rc) - np.array(get_row_column(m))))
            if dist < min_dist:
                min_dist = dist
                tbmt = m  # to be mutated to

        # itbm: index to be mutated
        itbm = random.sample([a.index for a in indi if a.symbol == tbmf[0]],
                             tbmf[1])
        for i in itbm:
            indi[i].symbol = tbmt
        syms = self.get_symbols_to_use(indi.get_chemical_symbols())
        indi.set_chemical_symbols(syms)
        if self.same_layer_composition:
            same_layer_comp(indi)

        parent_message = ': Parent {0}'.format(f.info['confid'])
        return (self.finalize_individual(indi),
                self.descriptor + parent_message)


class SymmetrySlabPermutation(SlabOperator):
    def __init__(self, verbose=False, num_muts=1, sym_goal=100, max_tries=50,
                 allowed_compositions=None, same_layer_composition=True):
        SlabOperator.__init__(self, verbose, num_muts,
                              allowed_compositions, same_layer_composition)
        try:
            import spglib
        except ImportError:
            print("SymmetrySlabPermutation needs spglib to function")

        assert sym_goal >= 1
        self.sym_goal = sym_goal
        self.max_tries = max_tries
        self.descriptor = 'SymmetrySlabPermutation'

    def get_new_individual(self, parents):
        f = parents[0]
        # Permutation only makes sense if two different elements are present
        if len(set(f.get_chemical_symbols())) == 1:
            f = parents[1]
            if len(set(f.get_chemical_symbols())) == 1:
                return None

        indi = self.initialize_individual(f, f)
        indi.info['data']['parents'] = [i.info['confid'] for i in parents]

        # Do the operation
        sym_num = 1
        sg = self.sym_goal
        while sym_num < sg:
            for _ in range(self.max_tries):
                for _ in range(2):
                    permute2(indi)
                if self.same_layer_composition:
                    same_layer_comp(indi)
                sym_num = spglib.get_symmetry_dataset(indi)['number']
                if sym_num >= sg:
                    break
            sg -= 1

        parent_message = ': Parent {0}'.format(f.info['confid'])
        return (self.finalize_individual(indi),
                self.descriptor + parent_message)


class RandomSlabPermutation(SlabOperator):
    def __init__(self, verbose=False, num_muts=1,
                 allowed_compositions=None, same_layer_composition=True):
        SlabOperator.__init__(self, verbose, num_muts,
                              allowed_compositions, same_layer_composition)

        self.descriptor = 'RandomSlabPermutation'

    def get_new_individual(self, parents):
        f = parents[0]
        # Permutation only makes sense if two different elements are present
        if len(set(f.get_chemical_symbols())) == 1:
            f = parents[1]
            if len(set(f.get_chemical_symbols())) == 1:
                return None

        indi = self.initialize_individual(f, f)
        indi.info['data']['parents'] = [i.info['confid'] for i in parents]

        # Do the operation
        for _ in range(self.num_muts):
            permute2(indi)
        if self.same_layer_composition:
            same_layer_comp(indi)

        parent_message = ': Parent {0}'.format(f.info['confid'])
        return (self.finalize_individual(indi),
                self.descriptor + parent_message)
