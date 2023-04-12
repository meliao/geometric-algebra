# Geometric-Algebra



## Overview of data representation

The object represents a multivector by storing its coefficients in the orthonormal basis
1, e_1, ..., e_n, e_1 e_2, ..., pseudoscalar. The coefficients are stored in an one-dimensional
array, and the ordering of this array is a bit counterintuitive. The ordering was chosen to simplify
the geometric product (hopefully).

Each basis element can be associated with a subset I of [n]. This subset can be interpreted as a 
bitstring of length n, where there are 1's in the indicies identified by I and 0's everywhere else. This bitstring
can be interpreted as a little-endian unsigned integer, which indicates the index of the basis element in the array.
```
index = bitstring[0] * (2 ** 0) + bitstring[1] * (2 ** 1) + bitstring[3] * (2 ** 3) + ... + bitstring[n] * (2 ** n)
```
Examples: 
 - e_1 corresponds to bitstring [1 0 0 ... 0] -> index 1
 - e_2 corresponds to bitstring [0 1 0 0 ... 0] -> index 2
 - e_3 corresponds to bitstring [0 0 1 0 ... 0] -> index 4
 - e_1 e_3 e_4 corresponds to bitstrinc [1 0 1 1 0 0 ...] -> index 1 + 4 + 8 = 13

## TODO

 <!-- - Rewrite `_get_zero_basis()` to return a flat 2^n length array -->
 <!-- - Rewrite `from_vector()` to put coefficients in positions (1, 2, 4, 8, ...) -->
 <!-- - Rewrite `add()` to operate on new data representation. -->
 <!-- - Re-make an ONB basis element product function which operates on integers, rather than binary strings. -->
 <!-- - Implement `geometric_product()` function -->
 - Implement `get_part(k)` function
 - Make a LUT for `_compute_output_index()` and `_compute_output_parity()` calls.



## Testing

```
python -m pytest test/
```
