/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	fp.h: Finite-field arithmetic in GF(p). All arithmetic in GF(p) is done
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

#ifndef FP_H
  #define FP_H

#include "fp2.h"
#include <gmp.h>
#include <math.h>

#include "sidh_params.h"

// Elements of GF(p^2)
typedef struct {
  mpz_t a;
} fp;

void fp_init(fp* x);

void fp_clear(fp *x);

int fp_setup(const char* characteristic);

void fp_free();

void fp_set_char(fp* x, const char* a);

void fp_set_mpz(fp* x, const mpz_t a);

void fp_set_characteristic(const char* p);

void fp_get(char *a, const fp x);

void fp_copy(fp* res, const fp x);

int fp_cmp(const fp x, const fp y);

int fp_is_one(const fp x);

int fp_is_zero(const fp x);

void fp_random(fp *res);

void fp_print(const fp x, char * a);

int fp_equals(fp *a, fp *b);

void fp_add(fp *res, const fp x, const fp y);

void fp_add_ui(fp *res, const fp x, unsigned long int u);

void fp_sub(fp *res, const fp x, const fp y);

void fp_sub_ui(fp *res, const fp x, unsigned long int u);

void fp_neg(fp *res, const fp x);

void fp_scalar(fp *res, const fp x, mpz_t s);

void fp_scalar_si(fp *res, const fp x, long int s);

void fp_mul(fp *res, const fp x, const fp y);

void fp_sqr(fp *res, const fp x);

int fp_inv(fp *res, const fp x);

int fp_div(fp *res, const fp x, const fp y);
  
#endif