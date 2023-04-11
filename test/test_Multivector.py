import numpy as np
import pytest

from src.Multivector import (MultiVector, ONBElement, _compute_num_transpositions,)


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

    
class Test_ONBElement:
    def test_0(self) -> None:
        """Easy case where the R side of the product just represents scalars in G_n
        """
        x = ONBElement(np.array([True, True, False]), 1)

        y = ONBElement(np.array([False, False, False]), 1)

        z = x.find_product(y)

        assert np.all(z.bool_array == np.array([True, True, False]))
        assert z.parity == 1

    def test_1(self) -> None:
        x = ONBElement(np.array([True, True, False]), 1)

        y = ONBElement(np.array([False, False, False]), -1)
        z = x.find_product(y)

        assert np.all(z.bool_array == np.array([True, True, False]))
        assert z.parity == -1

    def test_2(self) -> None:
        x = ONBElement(np.array([True, True, False, True]), 1)

        y = ONBElement(np.array([True, False, False, False]), 1)
        z = x.find_product(y)

        assert np.all(z.bool_array == np.array([False, True, False, True]))
        assert z.parity == 1

class Test_MultiVector:
    def test_0(self) -> None:
        a = np.array([2, 3, 4.5])
        b = np.array([1, 1, 1])

        mv_a = MultiVector.from_vector(a)

        assert len(mv_a.basis_expansion) == 4
        
        mv_b = MultiVector.from_vector(b)

        assert len(mv_b.basis_expansion) == 4

        mv_c = mv_a.add(mv_b)

        assert len(mv_c.basis_expansion) == 4
        assert np.allclose(mv_c.basis_expansion[1], a + b)