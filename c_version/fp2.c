/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	fp2.c: Finite-field arithmetic in GF(p^2). All arithmetic in GF(p) is done
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

#include "fp2.h"

void fp2_init(fp2* x) {
  mpz_init(x->a);
  mpz_init(x->b);
}

void fp2_clear(fp2 *x) {
  mpz_clear(x->a);
  mpz_clear(x->b);
}

// Initialization of GF(p,2)
int fp2_setup(const char* characteristic) {
  mpz_init( params_global.p );
  mpz_set_str(params_global.p, characteristic, 0);
  // Check that the Legendre symbol of -1 is -1 (p = 3 mod 4)
  if (mpz_fdiv_ui(params_global.p, 4) != 3) {
    mpz_clear(params_global.p);
    return 0;
  }

  gmp_randinit_default(params_global.state);
  mpz_init( params_global.randMult);
  rand_range(params_global.randMult,params_global.p);
  mpz_init(params_global.tmp1); mpz_init(params_global.tmp2); mpz_init(params_global.tmp3);
  int i;
  for (i = 0 ; i < GF_TMP_REGS ; i++)
    fp2_init(&params_global.fp2tmp[i]);
  fp2_init(&params_global.eqReg);
  fp2_init(&params_global.invReg);
  
  mpz_init(params_global.four_inverse);
  mpz_set_ui(params_global.four_inverse, 4);
  mpz_invert(params_global.four_inverse, params_global.four_inverse, params_global.p);
  
  params_global.initialized = 1;
  return 1;
}

void fp2_free(GF_params* field) {
  if (params_global.initialized) {
    int i;
    for (i = 0 ; i < GF_TMP_REGS ; i++)
      fp2_clear(&params_global.fp2tmp[i]);
    fp2_clear(&params_global.eqReg);
	fp2_clear(&params_global.invReg);
	mpz_clear(params_global.p); mpz_clear(params_global.tmp1); mpz_clear(params_global.tmp2); mpz_clear(params_global.tmp3);
	mpz_clear(params_global.randMult);
	mpz_clear(params_global.four_inverse);
	
    params_global.initialized = 0;
  }
}

// IO of elements
// outputs are strings in base 16
// (to save some bytes)
void fp2_set(fp2* x, const char* a, const char* b) {
  mpz_set_str(x->a, a, 0);
  mpz_set_str(x->b, b, 0);
}

void fp2_set_mpz(fp2* x, const mpz_t a, const mpz_t b){
  mpz_set(x->a,a);
  mpz_set(x->b,b);
}

void fp2_set_characteristic(const char* p) {
  mpz_set_str(params_global.p, p, 0);	
}

void fp2_get(char *a, char *b, const fp2 x) {
  gmp_sprintf(a, "%#Zx", x.a);
  gmp_sprintf(b, "%#Zx", x.b);
}

// Arithmetic modulo X^2 + 1
void fp2_copy(fp2* res, const fp2 x){
  mpz_set(res->a, x.a);
  mpz_set(res->b, x.b);
}

int fp2_cmp(const fp2 x, const fp2 y) {
  int c = mpz_cmp(x.a, y.a);
  if (c == 0) c = mpz_cmp(x.b, y.b);
  return c;
}

int fp2_is_one(const fp2 x) {
  return (mpz_sgn(x.a) == 0) && (mpz_cmp_ui(x.b, 1) == 0);
}

int fp2_is_zero(const fp2 x) {
  return (mpz_sgn(x.a) == 0) && (mpz_sgn(x.b) == 0);
}

void fp2_random(fp2 *res) {
  mpz_urandomm(res->a, params_global.state, params_global.p);
  mpz_urandomm(res->b, params_global.state, params_global.p);
}

void fp2_print(const fp2 x, char * a) {
  gmp_printf("%s: %Zd*x + %Zd\n",a, x.a, x.b);
}

//Returns 1 if a=b and another int otherwise
int fp2_equals(fp2 *a, fp2 *b){
  fp2_sub(&params_global.eqReg, *a, *b);
  return fp2_is_zero(params_global.eqReg);
}

void fp2_add(fp2 *res, const fp2 x, const fp2 y) {
  mpz_add(params_global.tmp1, x.a, y.a);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
  mpz_add(params_global.tmp1, x.b, y.b);
  mpz_mod(res->b, params_global.tmp1, params_global.p);
}

void fp2_add_ui(fp2 *res, const fp2 x, unsigned long int u) {
  mpz_add_ui(params_global.tmp1, x.b, u);
  mpz_mod(res->b, params_global.tmp1, params_global.p);
  mpz_set(res->a, x.a);
}

void fp2_sub(fp2 *res, const fp2 x, const fp2 y) {
  mpz_sub(params_global.tmp1, x.a, y.a);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
  mpz_sub(params_global.tmp1, x.b, y.b);
  mpz_mod(res->b, params_global.tmp1, params_global.p);
}

void fp2_sub_ui(fp2 *res, const fp2 x, unsigned long int u) {
  mpz_sub_ui(params_global.tmp1, x.b, u);
  mpz_mod(res->b, params_global.tmp1, params_global.p);
  mpz_set(res->a, x.a);
}

void fp2_neg(fp2 *res, const fp2 x) {
  if (mpz_sgn(x.a) == 0)
    mpz_set(res->a, x.a);
  else
    mpz_sub(res->a, params_global.p, x.a); 
  if (mpz_sgn(x.b) == 0)
    mpz_set(res->b, x.b);
  else
    mpz_sub(res->b, params_global.p, x.b);
}

void fp2_scalar(fp2 *res, const fp2 x, mpz_t s) {
  mpz_mul(params_global.tmp1, x.a, s);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
  mpz_mul(params_global.tmp1, x.b, s);
  mpz_mod(res->b, params_global.tmp1, params_global.p);
}

void fp2_scalar_si(fp2 *res, const fp2 x, long int s) {
  mpz_mul_si(params_global.tmp1, x.a, s);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
  mpz_mul_si(params_global.tmp1, x.b, s);
  mpz_mod(res->b, params_global.tmp1, params_global.p);
}

void fp2_mul(fp2 *res, const fp2 x, const fp2 y) {		
  mpz_add(params_global.tmp1, x.a, x.b);
  mpz_sub(params_global.tmp2, y.b, y.a);
  mpz_mul(params_global.tmp3, params_global.tmp1, params_global.tmp2);
  mpz_mul(params_global.tmp1, x.a, y.b);
  mpz_mul(params_global.tmp2, y.a, x.b);
  mpz_sub(params_global.tmp3, params_global.tmp3, params_global.tmp1);
  mpz_add(params_global.tmp1, params_global.tmp1, params_global.tmp2); 
  mpz_mod(res->a, params_global.tmp1, params_global.p);
  mpz_add(params_global.tmp3, params_global.tmp3, params_global.tmp2);
  mpz_mod(res->b, params_global.tmp3, params_global.p);
}

void fp2_sqr(fp2 *res, const fp2 x) {
  mpz_mul(params_global.tmp1, x.a, x.b);
  mpz_add(params_global.tmp1, params_global.tmp1, params_global.tmp1);
  mpz_add(params_global.tmp2, x.b, x.a);
  mpz_sub(params_global.tmp3, x.b, x.a);
  mpz_mul(params_global.tmp2, params_global.tmp2, params_global.tmp3);
  mpz_mod(res->a, params_global.tmp1, params_global.p);
  mpz_mod(res->b, params_global.tmp2, params_global.p);
}

int fp2_inv(fp2 *res, const fp2 x) {
  mpz_mul(params_global.tmp1, x.a, x.a);
  mpz_addmul(params_global.tmp1, x.b, x.b);
  mpz_mod(params_global.tmp1,params_global.tmp1,params_global.p);
  mpz_mul(params_global.tmp1,params_global.tmp1,params_global.randMult);
  mpz_mod(params_global.tmp1,params_global.tmp1, params_global.p);
  if (!mpz_invert(params_global.tmp3, params_global.tmp1, params_global.p))
	return 0;

  mpz_mul(params_global.tmp3,params_global.tmp3,params_global.randMult);
  mpz_mod(params_global.tmp3,params_global.tmp3,params_global.p);
  mpz_mul(params_global.tmp1, x.b, params_global.tmp3);
  mpz_neg(params_global.tmp3, params_global.tmp3);
  mpz_mul(params_global.tmp2, x.a, params_global.tmp3);

  mpz_mod(res->a, params_global.tmp2, params_global.p);
  mpz_mod(res->b, params_global.tmp1, params_global.p);

  return 1;
}

int fp2_div(fp2 *res, const fp2 x, const fp2 y) {
 if (!fp2_inv(&params_global.invReg, y)) return 0;
 fp2_mul(res, x, params_global.invReg);
 return 1;
}
