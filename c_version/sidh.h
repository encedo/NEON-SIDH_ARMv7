/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	sidh.h: Supersingular Isogeny Diffie-Hellman rounds and full key exchange
	
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

#ifndef SIDH_H
  #define SIDH_H
  
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <openssl/rand.h>

#include "fp.h"
#include "fp2.h"
#include "mont_curve.h"
#include "sidh_params.h"

void sidh_round1_alice(MP *Pother2, MP *Qother2, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, MP Pother, MP Qother, MP QP, MP PQ, int e);

void sidh_round1_bob(MP *Pother2, MP *Qother2, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, MP Pother, MP Qother, MP QP, MP PQ, int e);

void sidh_round2_alice(MC *E, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, int e);

void sidh_round2_bob(MC *E, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, int e);

   // returns 0 if shared keys are different, 1 if they are the same
int sidh_full(double *timeKey, double *timeAliceA, double *timeBobA, double *timeAliceB, double* timeBobB, char * eA, char * eB, char * lA_str, char * lB_str, int *strA, int lenA, int *strB, int lenB, MP *PA, MP *QA, MP *PB, MP *QB);
#endif