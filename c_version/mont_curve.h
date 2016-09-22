/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	mont_curve.h: Montgomery curve and point declarations and methods
	
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

#ifndef MONT_CURVE_H
  #define MONT_CURVE_H
  
#include "fp.h"
#include "fp2.h"

/* Utility routine to compute (A+2)/4 */
void MC_a24(fp2* A24, const fp2 A);

//Montgomery Curve
typedef struct {
  fp2 A, B, A24;
} MC;

void MC_copy(MC *res, const MC curve);

void MC_init(MC* x);

void MC_clear(MC *x);

void MC_print(const MC x);

void MC_set_curve( MC *curve, const fp2 A, const fp2 B, const fp2 A24 );

void MC_j_invariant(fp2 *final, const MC curve);

//Montgomery Point
typedef struct {
  fp2 x,y,z;
  MC curve;
} MP;

void MP_init(MP* a);
    
void MP_copy_point(MP *res, const MP P);
	
void MP_copy( MP *res, MP P );
    
void MP_set_curve( MP *res, const MC curve );

void MP_set( MP *res, const fp2 x, const fp2 y, const fp2 z, const MC curve );

void MP_clear(MP *a);

void MP_print(const MP a, char * string);

//P and Q must be on the same curve
void MP_add(MP *res, const MP P, const MP Q);

//P and Q must be on the same curve
void MP_sub(MP *res, const MP P, const MP Q);

//Returns the negative of P
void MP_neg(MP *res, const MP P);

void MP_x_affine(MP *res, const MP P);

//Kummer Point
typedef struct {
  fp2 x,z;
} KP;

void KP_init(KP* a);

void KP_clear(KP *a);
    
void KP_copy(KP *res, const fp2 x, const fp2 z);
	
void KP_copy_point(KP *res, const KP P);

void KP_set( KP *res, const fp2 x, const fp2 z);

void KP_print(const KP a, char * string);

#endif