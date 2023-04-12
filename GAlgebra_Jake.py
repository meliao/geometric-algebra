# GAlgebra code for geometric product

import numpy as np
from more_itertools import powerset as ps
import numbers

class GN:

    def __init__(self, n, epsilon=1.e-7):
        self.n = n
        self.epsilon=epsilon
        self.canonical_basis = ['1'] + [''.join(r) for i, r in enumerate(ps(['e' + str(j) for j in range(1, self.n+1)])) if not i == 0]
        self.where_m = np.zeros((2**self.n, 2**self.n))
        for i, b1 in enumerate(self.canonical_basis):
            for j, b2 in enumerate(self.canonical_basis):
                br, s = GN.basis_mult(b1, b2)
                if br == '1':
                    self.where_m[i,j] = s*self.epsilon
                else:
                    self.where_m[i,j] = s*self.canonical_basis.index(br)

    def __call__(self):
        return self.n

    def __eq__(self, other):
        if isinstance(other, GN):
            return self.n == other.n
        elif isinstance(other, numbers.Number):
            return self.n == other
        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, GN):
            return self.n < other.n
        elif isinstance(other, numbers.Number):
            return self.n < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, GN):
            return self.n <= other.n
        elif isinstance(other, numbers.Number):
            return self.n <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, GN):
            return self.n > other.n
        elif isinstance(other, numbers.Number):
            return self.n > other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, GN):
            return self.n >= other.n
        elif isinstance(other, numbers.Number):
            return self.n >= other
        else:
            return NotImplemented


    def basis_mult(b1, b2):
        if b1 == '1':
            return b2, 1
        if b2 == '1':
            return b1, 1
        l1, l2 = [int(i) for i in b1[1::2]], [int(i) for i in b2[1::2]]
        lr = set(l1).symmetric_difference(l2)
        matches, swaps = 0, 0
        for ind, i in enumerate(l1):
            if i in l2:
                swaps += ((len(l1) - 1) - ind) + (l2.index(i) - matches)
                matches += 1
        if swaps%2 == 0:
            s = 1
        else:
            s = -1
        if len(lr) == 0:
            return '1', s
        else:
            return ''.join(['e' + str(i) for i in lr]), s

    def mv(self, values):
        return MultiVector(self, values)

    def multivector(self, values):
        return MultiVector(self, values)


class MultiVector:

    def __init__(self, gn, values):
        """Create a multivector in Gn. Values should be given as either a list or a dictionary
        with entries representing the multivector in the canonical basis.
        Ex: MultiVector(3, {'1': 2, 'e1e3': 3, 'e1e2e3': -5}) = MultiVector(4, [2,0,0,0,0,3,0,-5])"""
        self.gn = gn
        if self.gn >= 10:
            raise ValueError("n is too large.")
        self.values = np.zeros(2**self.gn())
        if isinstance(values, list) or isinstance(values, np.ndarray):
            if not len(values) == 2**self.gn():
                raise ValueError("Incorrect number of values provided in list.")
            self.values = np.array(values)
        else:
            for p, v in values.items():
                self.values[self.gn.canonical_basis.index(p)] = v
        
    def __add__(self, other):
        if isinstance(other, MultiVector):
            return self.sum(other)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            return MultiVector(self.gn, other*self.values)
        if isinstance(other, MultiVector):
            return self.prod(other)
        return NotImplemented

    def sum(mv1, mv2):
        if not mv1.gn == mv2.gn:
            raise ValueError("Attempting to add multivectors from different spaces.")
        return MultiVector(mv1.gn, mv1.values + mv2.values)

    def prod(mv1, mv2):
        if not mv1.gn == mv2.gn:
            raise ValueError("Attempting to add multivectors from different spaces.")
        nv = [np.sum(np.where(np.isclose(abs(mv1.gn.where_m),i,atol=mv1.gn.epsilon*10),np.sign(mv1.gn.where_m)*np.outer(mv1.values,mv2.values),0)) for i in range(2**mv1.gn())]
        return MultiVector(mv1.gn, nv)

    def __str__(self):
        return str(self.values)

    def __repr__(self):
        return str(self.values)