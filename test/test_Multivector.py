import numpy as np
import pytest

from src.Multivector import (MultiVector, 
                            _compute_num_transpositions, 
                            _compute_output_index, 
                            _compute_output_parity,
                            _get_zero_basis)


class Test__compute_num_transpositions:
    def test_0(self) -> None:
        """
        Just tests to make sure everything runs without error
        """
        arr_0 = np.array([0,2,4,6])
        arr_1 = np.array([1,2,3,4])

        out = _compute_num_transpositions(arr_0, arr_1)

    def test_1(self) -> None:
        arr_0 = np.array([1, 1, 1, 1])
        arr_1 = np.array([1, 1])
        out = _compute_num_transpositions(arr_0, arr_1)

        assert out == 0

    def test_2(self) -> None:
        arr_0 = np.array([1, 3, 4])
        arr_1 = np.array([2, 3])
        out = _compute_num_transpositions(arr_0, arr_1)

        assert out == 3

class Test__compute_output_parity:
    def test_0(self) -> None:
        idx_1 = 13
        idx_2 = 11

        assert _compute_output_parity(idx_1, idx_2) == 1
        assert _compute_output_parity(idx_2, idx_1) == -1



class Test__compute_output_index:
    def test_0(self) -> None:

        for idx_1 in range(10):
            idx_0 = 0
            expected_out = idx_1

            out = _compute_output_index(idx_0, idx_1)
            assert out == expected_out

    def test_1(self) -> None:

        idx_0 = 13 # 1 0 1 1
        idx_1 = 11 # 1 1 0 1
        expected_out = 6 # 0 1 1 0

        out = _compute_output_index(idx_0, idx_1)
        assert out == expected_out

class Test_MultiVector:
    def test_0(self) -> None:
        a = np.array([2, 3, 4.5])
        b = np.array([1, 1, 1])

        mv_a = MultiVector.from_vector(a)

        assert len(mv_a.basis_expansion) == 8
        
        mv_b = MultiVector.from_vector(b)

        assert len(mv_b.basis_expansion) == 8

        mv_c = mv_a.add(mv_b)

        assert len(mv_c.basis_expansion) == 8

        expected_out = np.array([0.,
                                    2 + 1,
                                    3 + 1,
                                    0,
                                    4.5 + 1,
                                    0,
                                    0,
                                    0])
        assert np.allclose(mv_c.basis_expansion, expected_out)

    def test_1(self) -> None:
        """
        Tests geometric product of axis-aligned G_n elements 
        """

        basis_a = _get_zero_basis(3)
        basis_b = _get_zero_basis(3)

        # e_1
        basis_a[1] = 1.


        # e_3
        basis_b[4] = 1.

        mv_a = MultiVector(basis_a, 3)
        mv_b = MultiVector(basis_b, 3)

        # e_1 e_3
        mv_c = mv_a.geometric_product(mv_b)

        expected_basis = _get_zero_basis(3)
        expected_basis[5] = 1.

        assert np.all(mv_c.basis_expansion == expected_basis)

    def test_2(self) -> None:
        """
        Tests geometric product of two vectors in R_n
        """
        a = np.array([2, 3, 4.5])
        b = np.array([1, 1, 1])

        mv_a = MultiVector.from_vector(a)
        mv_b = MultiVector.from_vector(b)

        expected_basis = np.array([9.5, 0, 0, -1, 0, -2.5, -1.5, 0])

        mv_d = MultiVector(expected_basis, 3)

        mv_c = mv_a.geometric_product(mv_b)

        assert np.all(expected_basis == mv_c.basis_expansion)

        assert mv_d == mv_c


