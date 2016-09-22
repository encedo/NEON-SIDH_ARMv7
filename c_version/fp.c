/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	fp.c: Finite-field arithmetic in GF(p). All arithmetic in GF(p) is done
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

#include "fp.h"

GF_params params_global;

void fp_init(fp* x) {
  mpz_init(x->a);
}

void fp_clear(fp *x) {
  mpz_clear(x->a);
}

// IO of elements
// outputs are strings in base 16
// (to save some bytes)
void fp_set(fp* x, const char* a) {
  mpz_set_str(x->a, a, 0);
}

void fp_set_mpz(fp* x, const mpz_t a){
  mpz_set(x->a,a);
}

void fp_set_characteristic(const char* p) {
  mpz_set_str(params_global.p, p, 0);	
}

void fp_get(char *a, const fp x) {
  gmp_sprintf(a, "%#Zx", x.a);
}

// Arithmetic modulo X^2 + 1
void fp_copy(fp* res, const fp x){
  mpz_set(res->a, x.a);
}

int fp_cmp(const fp x, const fp y) {
  return mpz_cmp(x.a, y.a);
}

int fp_is_one(const fp x) {
  return (mpz_sgn(x.a) == 0);
}

int fp_is_zero(const fp x) {
  return (mpz_sgn(x.a) == 0);
}

void fp_random(fp *res) {
  mpz_urandomm(res->a, params_global.state, params_global.p);
}

void fp_print(const fp x, char * a) {
  gmp_printf("%s: %Zd\n",a, x.a);
}

//Returns 1 if a=b and another int otherwise
int fp_equals(fp *a, fp *b){
  fp_sub(&params_global.tmp1, *a, *b);
  return (mpz_sgn(params_global.tmp1) == 0);
}

void fp_add(fp *res, const fp x, const fp y) {
  mpz_add(params_global.tmp1, x.a, y.a);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_add_ui(fp *res, const fp x, unsigned long int u) {
  mpz_add_ui(params_global.tmp1, x.a, u);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_sub(fp *res, const fp x, const fp y) {
  mpz_sub(params_global.tmp1, x.a, y.a);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_sub_ui(fp *res, const fp x, unsigned long int u) {
  mpz_sub_ui(params_global.tmp1, x.a, u);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_neg(fp *res, const fp x) {
  if (mpz_sgn(x.a) == 0)
    mpz_set(res->a, x.a);
  else
    mpz_sub(res->a, params_global.p, x.a); 
}

void fp_scalar(fp *res, const fp x, mpz_t s) {
  mpz_mul(params_global.tmp1, x.a, s);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_scalar_si(fp *res, const fp x, long int s) {
  mpz_mul_si(params_global.tmp1, x.a, s);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_mul(fp *res, const fp x, const fp y) {		
  mpz_add(params_global.tmp1, x.a, y.a);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

void fp_sqr(fp *res, const fp x) {
  mpz_add(params_global.tmp1, x.a, x.a);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
}

int fp_inv(fp *res, const fp x) {
  if (!mpz_invert(res->a, x.a, params_global.p))
	return 0;
  return 1;
}

int fp_div(fp *res, const fp x, const fp y) {
  if (!mpz_invert(params_global.tmp1, y.a, params_global.p))
	return 0;
  mpz_mul(params_global.tmp2,params_global.tmp1,x.a);
  mpz_mod(res->a,params_global.tmp2,params_global.p);
  return 1;
}
