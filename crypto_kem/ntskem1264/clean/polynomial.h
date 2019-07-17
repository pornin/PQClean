/**
 *  polynomial.h
 *  NTS-KEM
 *
 *  Parameter: NTS-KEM(12, 64)
 *
 *  This file is part of the reference implemention of NTS-KEM
 *  submitted as part of NIST Post-Quantum Cryptography
 *  Standardization Process.
 **/

#ifndef _POLYNOMIAL_H
#define _POLYNOMIAL_H

#include <stdio.h>
#include <stdint.h>
#include "ff.h"

typedef struct {
    int size, degree;
    ff_unit *coeff;
} poly;

/**
 *  Initialise a new poly object
 *
 *  @param[in] size     The length of the polynomial
 *  @return The pointer to the polynomial object
 **/
poly* init_poly(int size);

/**
 *  Release a polynomial object
 *
 *  @param[in] px   The pointer to a polynomial object
 **/
void free_poly(poly* px);

/**
 *  Zero a polynomial
 *
 *  @param[in] px   The pointer to a polynomial object
 **/
void zero_poly(poly* px);

/**
 *  Clone a polynomial
 *
 *  @param[in] px   The pointer to a polynomial object
 *  @return The pointer to the cloned polynomial object
 **/
poly* clone_poly(const poly *px);

/**
 *  Compute the degree of a polynomial
 *
 *  @param[in] px   The pointer to a polynomial object
 **/
void update_poly_degree(poly *px);

/**
 *  Compute the formal derivative of a polynomial
 *
 *  @param[in]  fx  The pointer to a polynomial object
 *  @param[out] dx  The derivative of f(x)
 *  @return 1 on success, 0 otherwise
 **/
int formal_derivative_poly(const poly* fx, poly *dx);

/**
 *  Obtain the GCD of two finite-field polynomials
 *
 *  The method performs the following operation:
 *      g(x) = GCD(a(x), b(x))
 *
 *  @param ff2m The pointer to FF2m instance
 *  @param ax   The left input polynomial
 *  @param bx   The right input polynomial
 *  @param gx   The output polynomial
 *  @return 1 on successful operation, 0 otherwise
 **/
int gcd_poly(const FF2m* ff2m, const poly* ax, const poly *bx, poly *gx);

#if defined(INTERMEDIATE_VALUES)
/**
 *  Create a string representation of a polynomial
 *
 *  @param ff2m The pointer to FF2m instance
 *  @param ax   The input polynomial
 *  @return string representation on success, NULL otherwise
 **/
char* string_from_poly(const FF2m *ff2m, const poly *ax);

/**
 *  Print a polynomial to a file pointer
 *
 *  @param fp   The file pointer
 *  @param prefix   The prefix (if any)
 *  @param ff2m     The pointer to FF2m instance
 *  @param ax       The input polynomial
 *  @param suffix   The suffix (if any)
 **/
void poly_fprintf(FILE* fp,
                  const char *prefix,
                  const FF2m* ff2m,
                  const poly* ax,
                  const char *suffix);
#endif /* INTERMEDIATE_VALUES */

#endif /* _POLYNOMIAL_H */
