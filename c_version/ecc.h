/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	ecc.h: Elliptic curve cryptography (ECC) functions, primarily
		focused on Montgomery curves
	
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

#ifndef ECC_H
  #define ECC_H
  
#include "fp.h"
#include "fp2.h"
#include "mont_curve.h"

void mont_ladder(KP *R1, KP *R2, const KP P, const KP Q, const KP D, const fp2 A24);

void mont_double(KP *Q, const KP P, const fp2 A24);
		 
void mont_triple(KP *Q, const KP P, const fp2 A24);		 

// Three-point ladder addition step:
//   P1 = 2 P1
//   P2 = dadd(P1, P2, D2)
//   P3 = dadd(P1, P3, D3)		 
void mont_tradd(KP *P1, KP *P2, KP *P3, const fp2 dx2, const fp2 dx3, const fp2 A24);
		
void mont_3ladder(MP * R, const mpz_t t, const fp2 Px, const fp2 Qx, const fp2 QPx, const fp2 A24);
				  
void mont_to_ed(fp2* Rx, fp2* Ry, const MP P);
		
void shamir(MP *R, const MP P, const MP Q, const mpz_t m, const mpz_t n);		

#endif