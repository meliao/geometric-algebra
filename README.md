# Geometric-Algebra

## TODO

Design basis representation to align with binary strings. The XOR between two binary strings identifies the basis element created by the product of two G_n standard basis eleements. 

 <!-- - Rewrite `_get_zero_basis()` to return a flat 2^n length array -->
 - Rewrite `from_vector()` to put coefficients in positions (1, 2, 4, 8, ...)
 - Rewrite `add()` to operate on new data representation.
 - Re-make an ONB basis element product function which operates on integers, rather than binary strings.
 - Implement `geometric_product()` function
 - Implement `get_part(k)` function



## Testing

```
python -m pytest test/
```
