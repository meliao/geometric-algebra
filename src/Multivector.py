from tkinter import ON
from typing import List
import math

import numpy as np


class ONBElement:
    def __init__(self, bool_array: np.ndarray, parity: int) -> None:
        """An element in the orthonormal basis of G_n is represented as a boolean array
        which indicates the subset of R_n standard basis vectors that are producted together.
        I've also included a parity, because these products are only unique up to their sign, so
        I can represent both e_1 e_2 e_3 and e_2 e_1 e_3 with this data structure.  

        Args:
            bool_array (np.ndarray): Indicates the subset of R_n basis vectors in the G_n basis element 
                                      being represented.
            parity (int): Indicates the parity. This is either 1 or -1.
        """
        self.bool_array = bool_array
        self.parity = parity
        self.n = self.bool_array.shape[0]

    def find_product(self, other: None) -> None:
        """Given another ONBElement object, form their product: self * other
        Args:
            other (ONBElement): other ONB element to make the product
        """
        assert self.n == other.n

        # First, we find the xor. This determines the set of R_n basis vectors that'll be present 
        # in the product after canceling
        out_arr = np.logical_xor(self.bool_array, other.bool_array)

        # Second, number of transpositions to sort the concatenated arrays. This determines the 
        # parity of the product.

        # These two lines transform the boolean array into a list of integers which give the locations of the Trues 
        # in the boolean arrays
        idxes = np.nonzero(self.bool_array)[0]
        other_idxes = np.nonzero(other.bool_array)[0]

        n_transpositions = _compute_num_transpositions(idxes, other_idxes)
        if n_transpositions % 2:
            transposition_parity = -1
        else:
            transposition_parity = 1

        out_parity = self.parity * other.parity * transposition_parity
        return ONBElement(out_arr, out_parity)
        



def _compute_num_transpositions(arr_1: np.ndarray, arr_2: np.ndarray) -> int:
    """Computes the number of transpositions needed to sort arr_1 concatenated with arr_2.
    This might not be the optimal number of transpositions, but the parity will be correct.
    The parity of this number determines the parity of the product between two ONB elements.

    This function assumes arr_1 and arr_2 are sorted. It proceeds by determining the number of transpositions needed to 
    correctly insert arr_2[0] into arr_1, and then recursing.

    Args:
        arr_1 (np.ndarray): Sorted integer array
        arr_2 (np.ndarray): Sorted integer array

    Returns:
        int: Number of transpositions needed to sort the concatenation of arr_1 and arr_2
    """
    # Edge case where arr_2 is empty
    if arr_2.size == 0:
        return 0
    
    idx = np.searchsorted(arr_1, arr_2[0], side='right')
    num_transpositions = arr_1.shape[0] - idx

    if arr_2.shape[0] != 1:
        # Recursive case.

        return num_transpositions + _compute_num_transpositions(arr_1, arr_2[1:])

    else:
        # Base case when there's nothing in the arr_2 array
        return num_transpositions




def _get_zero_basis(k: int) -> List[np.ndarray]:
    """Given a max order k, this generates a list of arrays (with all zeros)
        that have the correct shape for an ONB of a k-vector.

    Args:
        k (int): max order

    Returns:
        List[np.ndarray]: List has length k + 1. Element j of this list is an array with length (k choose j)
    """
    out = []
    for i in range(k + 1):
        out.append(np.zeros(math.comb(k, i)))
    return out


class MultiVector:
    def __init__(self, basis_expansion: List[np.ndarray], n: int) -> None:
        self.basis_expansion = basis_expansion

        self.n = n

    @classmethod
    def from_blade(cls, blade: np.ndarray) -> None:
        raise NotImplementedError

    @classmethod
    def from_vector(cls, vector: np.ndarray) -> None:
        """Initializes a multivector in G_n from a vector in R_n

        Args:
            vector (np.ndarray): Has shape (n,)
        """
        empty_basis = _get_zero_basis(vector.shape[0])
        empty_basis[1] = vector
        return cls(empty_basis)


    def add(self, other: None) -> None:
        """Given another multivector, returns the sum of the two multivectors.

        self + other

        Args:
            other (MultiVector): The RHS of the sum
        """
        assert self.n == other.n


        out_basis = _get_zero_basis(self.n)

        for i in range(self.n):
            out_basis[i] += self.basis_expansion[i] + other.basis_expansion[i]
        
        return MultiVector(out_basis)

    def geometric_product(self, other: None) -> None:
        """Given another multivector, returns the geometric product of the 
        two multivectors.

        self * other

        Args:
            other (MultiVector): The RHS of the product
        """
        pass