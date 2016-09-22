/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	fp2.h: Finite-field arithmetic in GF(p^2). All arithmetic in GF(p) is done
		in the GNU Multiprecision Library (GMP)
	
	Copyright (c) Brian Koziel 2016

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FP2_H
  #define FP2_H

#include <gmp.h>
#include <math.h>

#include "sidh_params.h"

struct GF_params;
typedef struct GF_params GF_params;

// Elements of GF(p^2)
typedef struct {
  mpz_t a, b;
} fp2;

// The field GF(p^2)
// basically its characteristic and some work registers
struct GF_params {
  mpz_t p, tmp1, tmp2, tmp3, randMult, four_inverse;
  fp2 fp2tmp[GF_TMP_REGS];
  fp2 eqReg, invReg;
  gmp_randstate_t state;
  int initialized;
};

//Singleton pattern
extern GF_params params_global;

void fp2_init(fp2* x);

void fp2_clear(fp2 *x);

int fp2_setup(const char* characteristic);

void fp2_free();

void fp2_set_char(fp2* x, const char* a, const char* b);

void fp2_set_mpz(fp2* x, const mpz_t a, const mpz_t b);

void fp2_set_characteristic(const char* p);

void fp2_get(char *a, char *b, const fp2 x);

void fp2_copy(fp2* res, const fp2 x);

int fp2_cmp(const fp2 x, const fp2 y);

int fp2_is_one(const fp2 x);

int fp2_is_zero(const fp2 x);

void fp2_random(fp2 *res);

void fp2_print(const fp2 x, char * a);

int fp2_equals(fp2 *a, fp2 *b);

void fp2_add(fp2 *res, const fp2 x, const fp2 y);

void fp2_add_ui(fp2 *res, const fp2 x, unsigned long int u);

void fp2_sub(fp2 *res, const fp2 x, const fp2 y);

void fp2_sub_ui(fp2 *res, const fp2 x, unsigned long int u);

void fp2_neg(fp2 *res, const fp2 x);

void fp2_scalar(fp2 *res, const fp2 x, mpz_t s);

void fp2_scalar_si(fp2 *res, const fp2 x, long int s);

void fp2_mul(fp2 *res, const fp2 x, const fp2 y);

void fp2_sqr(fp2 *res, const fp2 x);

int fp2_inv(fp2 *res, const fp2 x);

int fp2_div(fp2 *res, const fp2 x, const fp2 y);

#endif