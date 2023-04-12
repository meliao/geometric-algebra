import numpy as np

        


def _compute_output_index(idx_1: int, idx_2: int) -> int:
    """Given two indices in the MultiVector's representation of the G_n orthonormal basis,
    this function computes the index of their output. 

    This function interprets the integers as little-endian bitstrings, takes their XOR, and
    interprets the result as an integer in little-endian.

    This function was written by chatgpt

    Args:
        idx_1 (int): input index 1
        idx_2 (int): input index 2

    Returns:
        int: output index
    """
    # bitwise XOR
    result = idx_1 ^ idx_2

    # Convert result to little endian format
    result_bytes = result.to_bytes((result.bit_length() + 7) // 8, byteorder='little')

    # Convert little endian bytes to integer
    result_int = int.from_bytes(result_bytes, byteorder='little')

    return result_int

def _compute_output_parity(idx_1: int, idx_2: int) -> int:
    """Given two indices in the MultiVector's representation of the G_n orthonormal basis,
    this function computes the parity of their output. The product of orthonormal basis elements
    in G_n result in (+1 / -1) of another basis element, and this function decides the parity of the product.

    The output parity is determined by taking the two integers, converting them to little endian bitstrings, 
    and then making arrays of their nonzero indices, and finding how many transpositions are required to 
    sort the concatenation of these arrays.

    Args:
        idx_1 (int): _description_
        idx_2 (int): _description_

    Returns:
        int: _description_
    """

    bitstring_1 = bin(idx_1)[2:][::-1]
    bitstring_2 = bin(idx_2)[2:][::-1]

    # Make sorted arrays containing the nonzero indices of the bitstrings
    idxes_1 = np.array([i for i, bit in enumerate(bitstring_1) if bit == '1'])
    idxes_2 = np.array([i for i, bit in enumerate(bitstring_2) if bit == '1'])

    # Count the number of transpositions required to sort the concatenation
    # [idxes_1, idxes_2]

    n_transpositions = _compute_num_transpositions(idxes_1, idxes_2)

    if n_transpositions % 2:
        return -1

    else:
        return 1





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




def _get_zero_basis(k: int) -> np.ndarray:
    """Lterally just initializes an array filled with zeros of length 2**k

    Args:
        k (int): max order

    Returns:
       np.ndarray: Has length 2**k
    """
    return np.zeros(2**k)


class MultiVector:
    def __init__(self, basis_expansion: np.array, n: int) -> None:
        """Represents a multivector in G_n, and implements addition and geometric product. 

        The object represents a multivector by storing its coefficients in the orthonormal basis
        1, e_1, ..., e_n, e_1 e_2, ..., pseudoscalar. The coefficients are stored in an one-dimensional
        array, and the ordering of this array is a bit counterintuitive. The ordering was chosen to simplify
        the geometric product (hopefully).

        Each basis element can be associated with a subset I of [n]. This subset can be interpreted as a 
        bitstring of length n, where there are 1's in the indicies identified by I and 0's everywhere else. This bitstring
        can be interpreted as a little-endian unsigned integer, which indicates the index of the basis element in the array.
        
        index = bitstring[0] * (2 ** 0) + bitstring[1] * (2 ** 1) + bitstring[3] * (2 ** 3) + ... + bitstring[n] * (2 ** n)

        So e_1 corresponds to bitstring [1 0 0 ... 0] -> index 1

        e_2 corresponds to bitstring [0 1 0 0 ... 0] -> index 2

        e_3 corresponds to bitstring [0 0 1 0 ... 0] -> index 4

        e_1 e_3 e_4 corresponds to bitstrinc [1 0 1 1 0 0 ...] -> index 1 + 4 + 8 = 13

        Args:
            basis_expansion (np.array): Basis expansion with the ordering mentioned above. Has shape (2**n)
            n (int): Corresponds to the n in G_n
        """
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
        n = vector.shape[0]
        empty_basis = _get_zero_basis(n)

        for i in range(n):
            empty_basis[2**i] = vector[i]
        
        return cls(empty_basis, n)


    def add(self, other: None) -> None:
        """Given another multivector, returns the sum of the two multivectors.

        self + other

        Args:
            other (MultiVector): The RHS of the sum
        """
        assert self.n == other.n

        return MultiVector(self.basis_expansion + other.basis_expansion, self.n)


    def geometric_product(self, other: None) -> None:
        """Given another multivector, returns the geometric product of the 
        two multivectors.

        self * other

        Args:
            other (MultiVector): The RHS of the product
        """
        assert self.n == other.n
        
        out_basis = _get_zero_basis(self.n)

        for i in range(2 ** self.n):
            for j in range(2 ** self.n):

                # This could be a LUT if I were motivated
                output_idx = _compute_output_index(i, j)
                output_parity = _compute_output_parity(i, j)

                out_basis[output_idx] += self.basis_expansion[i] * other.basis_expansion[j] * output_parity

        return MultiVector(out_basis, self.n)


    def __eq__(self, other) -> bool:
        """
        Evaluates equality between two multivectors
        """

        return self.n == other.n and np.all(self.basis_expansion == other.basis_expansion)

    def __ne__(self, other) -> bool:
        """
        Evaluates inequality between two multivectors
        """
        return not self.__eq__(other)