/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	ecc.c: Elliptic curve cryptography (ECC) functions, primarily
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

#include "ecc.h"
// One step of Montgomery ladder
void mont_ladder(KP *R1,
		 KP *R2,
		 const KP P,
		 const KP Q,
		 const KP D, //Differential
		 const fp2 A24) {
  fp2* tmp = params_global.fp2tmp;

  fp2_add(&tmp[4], P.x, P.z);         // a = (self.x + self.z)
  fp2_sub(&tmp[5], P.x, P.z);         // b = (self.x - self.z)
  fp2_sub(&tmp[6], Q.x, Q.z);
  fp2_add(&tmp[7], Q.x, Q.z);
  fp2_sqr(&tmp[1], tmp[4]);         // aa = a.square()
  fp2_sqr(&tmp[3], tmp[5]);         // bb = b.square()
  fp2_sub(&tmp[0], tmp[1], tmp[3]); // e = aa - bb
  fp2_mul(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
  fp2_mul(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b

  fp2_add(&tmp[2], tmp[6], tmp[7]);
  fp2_sqr(&tmp[2], tmp[2]);
  fp2_mul(&tmp[2], tmp[2], D.z);     // x2 = diff.z*(da + cb).square()

  fp2_sub(&tmp[8], tmp[6], tmp[7]);
  fp2_sqr(&tmp[8], tmp[8]);
  fp2_mul(&tmp[8], tmp[8], D.x);     // z2 = diff.x*(da - cb).square()

  fp2_mul(&(R1->z), A24, tmp[0]);
  fp2_add(&(R1->z), R1->z, tmp[3]);
  fp2_mul(&(R1->z), R1->z, tmp[0]);   // z1 = e*(bb + self.curve.A24*e))
  fp2_mul(&(R1->x), tmp[1], tmp[3]);  // x1 = aa*bb

  fp2_copy(&(R2->x), tmp[2]);
  fp2_copy(&(R2->z), tmp[8]);
}

/* Montgomery point doubling */
void mont_double(KP *Q,
		 const KP P,
		 const fp2 A24) {

  fp2* tmp = params_global.fp2tmp;
  fp2_add(&tmp[0], P.x, P.z);           // a = (x + z)
  fp2_sqr(&tmp[1], tmp[0]);         // aa = a^2
  fp2_sub(&tmp[2], P.x, P.z);           // b = (x - z)
  fp2_sqr(&tmp[3], tmp[2]);         // bb = b^2
  fp2_sub(&tmp[4], tmp[1], tmp[3]); // c = aa - bb
  fp2_mul(&(Q->z), A24, tmp[4]);   
  fp2_add(&(Q->z), Q->z, tmp[3]);			
  fp2_mul(&(Q->z), Q->z, tmp[4]);   // Z = c (bb + A24 c))	
  fp2_mul(&(Q->x), tmp[1], tmp[3]);  // X = aa bb  
}

/* Montgomery point tripling */
void mont_triple(KP * Q,
		 const KP P,
		 const fp2 A24) {
  fp2* tmp = params_global.fp2tmp;
  
  fp2_add(&tmp[0], P.x, P.z);           // a = (x + z)
  fp2_sqr(&tmp[1], tmp[0]);         // aa = a^2
  fp2_sub(&tmp[2], P.x, P.z);           // b = (x - z)
  fp2_sqr(&tmp[3], tmp[2]);         // bb = b^2
  fp2_sub(&tmp[4], tmp[1], tmp[3]); // c = aa - bb
  fp2_mul(&tmp[6], A24, tmp[4]);   
  fp2_add(&tmp[6], tmp[6], tmp[3]);			
  fp2_mul(&tmp[6], tmp[6], tmp[4]);   // Z = c (bb + A24 c))	
  fp2_mul(&tmp[5], tmp[1], tmp[3]);  // X = aa bb  

  fp2_sub(&tmp[7], tmp[5], tmp[6]);
  fp2_add(&tmp[8], tmp[5], tmp[6]);
  fp2_mul(&tmp[5], tmp[7], tmp[0]); // da = (x2 - z2)*a
  fp2_mul(&tmp[6], tmp[8], tmp[2]); // cb = (x2 + z2)*b

  fp2_add(&tmp[7], tmp[5], tmp[6]);
  fp2_sqr(&tmp[7], tmp[7]);
  fp2_mul(&tmp[7], tmp[7], P.z);     // X = z*(da + cb)^2

  fp2_sub(&tmp[8], tmp[5], tmp[6]);
  fp2_sqr(&tmp[8], tmp[8]);
  fp2_mul(&(Q->z), tmp[8], P.x);     // Z = x*(da - cb)^2

  //Copy to ensure P.x isn't overwritten
  fp2_copy(&(Q->x), tmp[7]);
}

// Three-point ladder addition step:
//   P1 = 2 P1
//   P2 = dadd(P1, P2, D2)
//   P3 = dadd(P1, P3, D3)
void mont_tradd(KP *P1,
        KP *P2,
		KP *P3,
		const fp2 dx2,
		const fp2 dx3,
		const fp2 A24) {
    fp2* tmp = params_global.fp2tmp;
	
    fp2_add(&tmp[4], P1->x, P1->z);         // a = (self.x + self.z)
    fp2_sub(&tmp[5], P1->x, P1->z);         // b = (self.x - self.z)
    fp2_sub(&tmp[6], P2->x, P2->z);
    fp2_add(&tmp[7], P2->x, P2->z);
    fp2_sub(&tmp[0], P3->x, P3->z);
    fp2_add(&tmp[1], P3->x, P3->z); 
    
    fp2_mul(&tmp[2], tmp[0], tmp[4]); // da = (P.x - P.z)*a
    fp2_mul(&tmp[3], tmp[1], tmp[5]); // cb = (P.x + P.z)*b
    

    fp2_add(&tmp[8], tmp[2], tmp[3]);
	fp2_sqr(&(P3->x), tmp[8]);
    
    fp2_sub(&tmp[8], tmp[2], tmp[3]);
    fp2_sqr(&tmp[0], tmp[8]);
    fp2_mul(&(P3->z), tmp[0], dx3);     // z2 = diff.x*(da - cb).square()
    
    /* P2 */
    fp2_mul(&tmp[6], tmp[6], tmp[4]); // da = (P.x - P.z)*a
    fp2_mul(&tmp[7], tmp[7], tmp[5]); // cb = (P.x + P.z)*b
    
    fp2_add(&tmp[2], tmp[6], tmp[7]);
    fp2_sqr(&(P2->x), tmp[2]);
    
    fp2_sub(&tmp[3], tmp[6], tmp[7]);
    fp2_sqr(&tmp[8], tmp[3]);
    fp2_mul(&(P2->z), tmp[8], dx2);     // z2 = diff.x*(da - cb).square()
    
    
    /* P1 */
    fp2_sqr(&tmp[6], tmp[4]);         // aa = a.square()
    fp2_sqr(&tmp[7], tmp[5]);         // bb = b.square()
    fp2_sub(&tmp[8], tmp[6], tmp[7]); // e = aa - bb
    fp2_mul(&tmp[4], A24, tmp[8]);
    fp2_add(&tmp[5], tmp[4], tmp[7]);
    fp2_mul(&(P1->z), tmp[5], tmp[8]);      // z1 = e*(bb + self.curve.A24*e))
    fp2_mul(&(P1->x), tmp[6], tmp[7]);      // x1 = aa*bb
   
}


// 3-point ladder to compute P + [t]Q
// Inputs: t, P, Q, Q - P
void mont_3ladder(MP *R,
                  const mpz_t t,
                  const fp2 Px,
                  const fp2 Qx,
                  const fp2 QPx,
                  const fp2 A24) {
  fp2* Rx = &(R->x); fp2* Rz = &(R->z);
  fp2* tmp = params_global.fp2tmp;
  KP *P1 = malloc(sizeof(KP));
  KP *P2 = malloc(sizeof(KP));
  KP *P3 = malloc(sizeof(KP));
  
  KP_init(P1); KP_init(P2); KP_init(P3);
  
  fp2_set(&(P1->x), "0", "1");
  fp2_set(&(P1->z), "0", "0");
  
  fp2_copy(&(P2->x), Qx);
  fp2_copy(&(P2->z), P1->x);
  
  fp2_copy(&(P3->x), Px); 
  fp2_copy(&(P3->z), P1->x);
  
  int bit = mpz_sizeinbase(t, 2) - 1;
  for ( ; bit >=0 ; bit--) {
    if (mpz_tstbit(t, bit) == 0) {
       mont_tradd(P1, P2, P3,
                  Qx, Px, A24);
    } else {
       mont_tradd(P2, P1, P3,
                  Qx, QPx, A24);
    }
  }
  
  fp2_copy(Rx,P3->x); fp2_copy(Rz,P3->z);
  KP_clear(P1); KP_clear(P2); KP_clear(P3);
  free(P1); free(P2); free(P3);
}

void mont_to_ed(fp2* Rx, fp2* Ry,
		const MP P) {
  fp2* tmp = params_global.fp2tmp;
  /*
    X = x(x+z) / y(x+z)
    Y = y(x-z) / y(x+z)
  */
  fp2_add(&tmp[5], P.x, P.z);
  fp2_sub(&tmp[6], P.x, P.z);
  fp2_mul(&tmp[7], P.y, tmp[5]);
  fp2_inv(&tmp[8], tmp[7]);
  fp2_mul(&tmp[7], tmp[5], tmp[8]);
  fp2_mul(Rx, P.x, tmp[7]);
  fp2_mul(&tmp[7], tmp[6], tmp[8]);
  fp2_mul(Ry, P.y, tmp[7]);
}
/*
  Computes [m]P + [n]Q, with P and Q points on the Montgomery curve
  with parameters A,B.  Uses Edwards' coordinates for
  calculations.  */
void shamir(MP *R,
	    const MP P,
	    const MP Q,
	    const mpz_t m, const mpz_t n) {
  // some temporary registers
  fp2* tmp = params_global.fp2tmp;
  fp2* Rx = &(R->x); fp2* Ry = &(R->y); fp2* Rz = &(R->z);
  
  fp2_inv(&tmp[0], P.curve.B);
  fp2_add_ui(&tmp[1], P.curve.A, 2);
  fp2_mul(&tmp[10], tmp[1], tmp[0]);
  fp2_sub_ui(&tmp[1], tmp[1], 4);
  fp2_mul(&tmp[11], tmp[1], tmp[0]);

  mont_to_ed(&tmp[12], &tmp[13], P);
  mont_to_ed(&tmp[14], &tmp[15], Q);

  /*
    Computing P+Q using affine Edwards.
  */
  fp2_mul(&tmp[4], tmp[12], tmp[14]); // tmp4 = C = aPx * aQx
  fp2_mul(&tmp[5], tmp[13], tmp[15]); // tmp5 = D = aPy * aQy
  fp2_add(&tmp[0], tmp[12], tmp[13]); // tmp0 = A = aPx + aPy
  fp2_add(&tmp[2], tmp[14], tmp[15]); // tmp2 = B = aQx + aQy
  fp2_mul(&tmp[7], tmp[4], tmp[5]);
  fp2_mul(&tmp[6], tmp[11], tmp[7]); // tmp6 = E = d * aPx * aQx * aPy * aQy
  fp2_sqr(&tmp[8], tmp[6]);
  fp2_neg(&tmp[7], tmp[8]);
  fp2_add_ui(&tmp[7], tmp[7], 1);
  fp2_inv(&tmp[8], tmp[7]); // tmp8 = 1 / (1-E^2)
  fp2_add_ui(&tmp[6], tmp[6], 1);
  fp2_mul(&tmp[7], tmp[10], tmp[4]);
  fp2_sub(&tmp[1], tmp[5], tmp[7]);
  fp2_mul(&tmp[7], tmp[6], tmp[1]);
  fp2_mul(&tmp[9], tmp[7], tmp[8]); // PQy = (1+E)(D - a C) / (1-E^2)
  fp2_neg(&tmp[6], tmp[6]);
  fp2_add_ui(&tmp[6], tmp[6], 2);
  fp2_mul(&tmp[1], tmp[0], tmp[2]);
  fp2_sub(&tmp[3], tmp[1], tmp[4]);
  fp2_sub(&tmp[1], tmp[3], tmp[5]);
  fp2_mul(&tmp[7], tmp[6], tmp[1]);
  fp2_mul(&tmp[16], tmp[7], tmp[8]); // PQx = (1-E)(A B - C - D) / (1-E^2)*/
  
  int bit = MAX(mpz_sizeinbase(m, 2), mpz_sizeinbase(n, 2)) - 1;
  mpz_set_ui(Rx->a, 0); mpz_set_ui(Ry->a, 0); mpz_set_ui(Rz->a, 0);   
  mpz_set_ui(Rx->b, 0); mpz_set_ui(Ry->b, 1); mpz_set_ui(Rz->b, 1);

  for ( ; bit >=0 ; bit--){
    /* Double, using projective Edwards */
    fp2_add(&tmp[1], *Rx, *Ry);
    fp2_sqr(&tmp[0], tmp[1]); // tmp0 = B = (Rx + Ry)^2
    fp2_sqr(&tmp[1], *Rx); // tmp1 = C = Rx^2
    fp2_sqr(&tmp[2], *Ry); // tmp2 = D = Ry^2
    fp2_mul(&tmp[3], tmp[10], tmp[1]); // tmp3 = E = a C
    fp2_add(&tmp[4], tmp[3], tmp[2]); // tmp4 = F = E + D
	fp2_sqr(&tmp[5], *Rz); // tmp5 = H = Rz^2
    fp2_scalar_si(&tmp[7], tmp[5], 2);
    fp2_sub(&tmp[6], tmp[4], tmp[7]); // tmp6 = J = F - 2H
    fp2_sub(&tmp[7], tmp[0], tmp[1]);
    fp2_sub(&tmp[8], tmp[7], tmp[2]);
    fp2_mul(Rx, tmp[8], tmp[6]); // Rx = (B-C-D) J
    fp2_sub(&tmp[7], tmp[3], tmp[2]);
    fp2_mul(Ry, tmp[7], tmp[4]); // Ry = (E-D) F
    fp2_mul(Rz, tmp[4], tmp[6]); // Rz = F J

    /* Double and Add, using projective Edwards */
    int r = mpz_tstbit(m, bit) | (mpz_tstbit(n, bit) << 1);
    if (r) {
      if (r == 1) {
	fp2_mul(&tmp[0], *Rx, tmp[12]); // tmp0 = C = Rx aPx
	fp2_mul(&tmp[1], *Ry, tmp[13]); // tmp1 = D = Ry aPy
	fp2_add(&tmp[2], tmp[12], tmp[13]);  // tmp2 = H = aPx + aPy
      } else if (r == 2) {
	fp2_mul(&tmp[0], *Rx, tmp[14]); // tmp0 = C = Rx aQx
	fp2_mul(&tmp[1], *Ry, tmp[15]); // tmp1 = D = Ry aQy
	fp2_add(&tmp[2], tmp[14], tmp[15]);  // tmp2 = H = aQx + aQy
      } else {
	fp2_mul(&tmp[0], *Rx, tmp[16]); // tmp0 = C = Rx PQx
	fp2_mul(&tmp[1], *Ry, tmp[9]); // tmp1 = D = Ry PQy
	fp2_add(&tmp[2], tmp[16], tmp[9]);  // tmp2 = H = PQx + PQy
      }
      fp2_sqr(&tmp[3], *Rz); // tmp3 = B = Rz^2
      fp2_mul(&tmp[5], tmp[0], tmp[1]);
      fp2_mul(&tmp[4], tmp[11], tmp[5]); // tmp4 = E = d C D
      fp2_sub(&tmp[5], tmp[3], tmp[4]); // tmp5 = F = B - E
      fp2_add(&tmp[6], tmp[3], tmp[4]); // tmp6 = G = B + E
      fp2_add(&tmp[7], *Rx, *Ry);
      fp2_mul(&tmp[8], tmp[7], tmp[2]);
      fp2_sub(&tmp[7], tmp[8], tmp[0]);
      fp2_sub(&tmp[8], tmp[7], tmp[1]);
      fp2_mul(&tmp[7], tmp[8], tmp[5]);
      fp2_mul(Rx, *Rz, tmp[7]); // Rx = Rz F ((Rx+Ry)H - C - D)
      fp2_mul(&tmp[7], tmp[10], tmp[0]);
      fp2_sub(&tmp[8], tmp[1], tmp[7]);
      fp2_mul(&tmp[7], tmp[6], tmp[8]);
      fp2_mul(Ry, *Rz, tmp[7]); // Ry = Rz G (D - a C)
      fp2_mul(Rz, tmp[5], tmp[6]); // Rz = F G
    }
  }
  /* Convert to Montgomery */
  fp2_add(&tmp[0], *Rz, *Ry);
  fp2_sub(&tmp[1], *Rz, *Ry);
  fp2_mul(&(R->y), tmp[0], *Rz); // Ry = (Rz+Ry)Rz
  fp2_mul(&(R->z), tmp[1], *Rx); // Rz = (Rz-Ry)Rx
  fp2_mul(&(R->x), tmp[0], *Rx); // Rx = (Rz+Ry)Rx

}
