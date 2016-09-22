/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	iso.c: Isogeny structures and computational methods
	
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

#include "iso.h"

/*
  Compute an isomorphism of the montgomery curve
  sending (x,z) to (0,0).
*/
void isom_comp(iso* iso, MC *phiE, 
           const MC E,
	       const KP P) {
  fp2* tmp = params_global.fp2tmp;
  
  //mont_double(&tmp[1], &tmp[2], P.x, P.z, A24);
  fp2_add(&tmp[0], P.x, P.z);           // a = (x + z)
  fp2_sqr(&tmp[1], tmp[0]);         // aa = a^2
  fp2_sub(&tmp[2], P.x, P.z);           // b = (x - z)
  fp2_sqr(&tmp[3], tmp[2]);         // bb = b^2
  fp2_sub(&tmp[4], tmp[1], tmp[3]); // c = aa - bb
  fp2_mul(&tmp[2], E.A24, tmp[4]);   
  fp2_add(&tmp[2], tmp[2], tmp[3]);			
  fp2_mul(&tmp[2], tmp[2], tmp[4]);   // Z = c (bb + A24 c))	
  fp2_mul(&tmp[1], tmp[1], tmp[3]);  // X = aa bb  
  
  fp2_div(&tmp[0], tmp[1], tmp[2]);  // P2x = x([2]P) / z([2]P)
  fp2_neg(&iso->r, tmp[0]);  // r = -P2x
  fp2_scalar_si(&tmp[1], tmp[0], 3);
  fp2_add(&tmp[1], tmp[1], E.A);  // a2 = 3 P2x + A
  fp2_mul(&tmp[2], iso->r, P.z);
  fp2_add(&tmp[3], tmp[2], P.x);
  fp2_div(&iso->u, P.z, tmp[3]);  // u = z / (x - z P2x)
  fp2_mul(&(phiE->A), tmp[1], iso->u);  // iA = a2 u
  fp2_mul(&(phiE->B), E.B, iso->u);  // iB = B u
  MC_a24(&(phiE->A24), phiE->A);
}

/* Apply an isomorphism of Montgomery curves */
void isom_eval_kummer(KP *Q,
		const iso iso,
		const KP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_mul(&tmp[0], iso.r, P.z);
  fp2_add(&tmp[1], P.x, tmp[0]);
  fp2_mul(&(Q->x), tmp[1], iso.u); // X = (x + r z) u
}

/* Apply an isomorphism of Montgomery curves */
void isom_eval_proj(MP *Q,
		const iso iso,
		const MP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_mul(&tmp[0], iso.r, P.z);
  fp2_add(&tmp[1], P.x, tmp[0]);
  fp2_mul(&(Q->y), P.y, iso.u); // Y = y u
  fp2_mul(&(Q->x), tmp[1], iso.u); // X = (x + r z) u
}

/*
  Compute a 2-isogeny of the montgomery curve
  sending (x,z) to (1,...).
*/
void iso2_comp(iso2* iso, MC * phiE,
	       const MC E,
	       const KP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_sub(&tmp[0], P.x, P.z);
  fp2_sqr(&tmp[1], tmp[0]);
  fp2_inv(&tmp[0], tmp[1]);
  fp2_mul(&tmp[1], tmp[0], P.z);
  fp2_mul(iso, tmp[1], P.x); // iA2 = x z / (x-z)^2
  fp2_add_ui(&tmp[0], E.A, 6);
  fp2_mul(&(phiE->B), E.B, *iso); // iB = B iA2
  fp2_mul(&(phiE->A), tmp[0], *iso); // iA = (A+6) iA2
  MC_a24(&(phiE->A24), phiE->A);
}

/* Apply a 2-isogeny of Montgomery curves */
void iso2_eval_kummer(KP *Q,
		const iso2 iso,
		const KP P) {
  fp2* tmp = params_global.fp2tmp;
  
  fp2_sub(&tmp[3], P.x, P.z);
  fp2_sqr(&tmp[4], tmp[3]);
  fp2_mul(&(Q->z), P.z, P.x); // Z = z x
  fp2_mul(&(Q->x), iso, tmp[4]); // X = iA2 (x - z)^2
}

/* Apply a 2-isogeny of Montgomery curves */
void iso2_eval_proj(MP *Q,
		const iso2 iso,
		const MP P) {
  fp2* tmp = params_global.fp2tmp;
  
  fp2_sub(&tmp[3], P.x, P.z);
  fp2_sqr(&tmp[4], tmp[3]);
  fp2_mul(&tmp[4], P.x, tmp[4]); // ... X = x iA2 (x - z)^2
  fp2_sqr(&tmp[0], P.x); // Px2 = x^2
  fp2_sqr(&tmp[1], P.z);
  fp2_sub(&tmp[2], tmp[0], tmp[1]);
  fp2_mul(&tmp[1], P.y, tmp[2]);
  fp2_mul(&(Q->y), iso, tmp[1]); // Y = iA2 y (x^2 - z^2)
  fp2_mul(&(Q->z), P.z, tmp[0]); // Z = z x^2
  fp2_mul(&(Q->x), iso, tmp[4]); // X = iA2 (x - z)^2
}

/*
  Compute a 3-isogeny of the montgomery curve
*/
void iso3_comp(iso3* iso, MC *phiE,
	       const MC E,
	       const KP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_div(&iso->p, P.x, P.z);             // p
  fp2_sqr(&iso->p2, iso->p);          // p^2
  
  fp2_scalar_si(&tmp[3], iso->p, -6); 
  fp2_add(&tmp[4], tmp[3], E.A);
  fp2_mul(&tmp[3], tmp[4], iso->p);
  fp2_add_ui(&tmp[4], tmp[3], 6);     // (-6p + A)p + 6

  fp2_mul(&(phiE->B), E.B, iso->p2);      // iB = B p^2
  fp2_mul(&(phiE->A), tmp[4], iso->p);  // iA = ((-6p + A)p + 6)p

  MC_a24(&(phiE->A24), phiE->A);
}

/* Apply a 3-isogeny of Montgomery curves */
void iso3_eval_kummer(KP *Q,
		const iso3 iso,
		const KP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_mul(&tmp[0], P.z, iso.p);
  fp2_sub(&tmp[1], P.x, tmp[0]); // h = x - p z
                              // if zero, P is in the kernel
  fp2_mul(&tmp[2], P.x, iso.p);
  fp2_sub(&tmp[0], tmp[2], P.z); // rh = x p - z
  fp2_sqr(&tmp[3], tmp[0]);
  fp2_mul(&tmp[2], P.x, tmp[3]); // X0 = x (x p - z)^2

  fp2_sqr(&tmp[3], tmp[1]);
  fp2_mul(&(Q->z), tmp[3], P.z); // Z = h^2 z
  fp2_copy(&(Q->x), tmp[2]);
}
/* Apply a 3-isogeny of Montgomery curves */
void iso3_eval_proj(MP *Q,
		const iso3 iso,
		const MP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_mul(&tmp[0], P.z, iso.p);
  fp2_sub(&tmp[1], P.x, tmp[0]); // h = x - p z
                              // if zero, P is in the kernel
  fp2_mul(&tmp[2], P.x, iso.p);
  fp2_sub(&tmp[0], tmp[2], P.z); // rh = x p - z
  fp2_sqr(&tmp[3], tmp[0]);
  fp2_mul(&tmp[2], P.x, tmp[3]); // X0 = x (x p - z)^2
  fp2_mul(&tmp[2], tmp[2], tmp[1]); // X0 *= h
  fp2_mul(&tmp[3], P.x, P.z);
  fp2_sub_ui(&tmp[4], iso.p2, 1);
  fp2_mul(&tmp[5], tmp[3], tmp[4]);
  fp2_scalar_si(&tmp[3], tmp[5], -2);
  fp2_mul(&tmp[4], tmp[0], tmp[1]);
  fp2_add(&tmp[5], tmp[3], tmp[4]);
  fp2_mul(&tmp[3], tmp[5], tmp[0]); // (rh (rh h + 2xz(1-p^2)))
  fp2_sqr(&tmp[7], tmp[1]);
  fp2_mul(&tmp[8], tmp[7], P.z);
  fp2_mul(&(Q->y), P.y, tmp[3]); // Y = y (rh (rh h + 2xz(1-p^2)))
  fp2_mul(&(Q->z), tmp[8], tmp[1]); // Z = h^2 h z
  fp2_copy(&(Q->x), tmp[2]);
}



/*
  Compute a 4-isogeny of the Montgomery curve
  sending (1,...) to infinity.
*/
void iso4_comp(iso4* iso, MC *phiE,
	       const MC E) {
  fp2* tmp = params_global.fp2tmp;

  fp2_add_ui(&iso->Ap2, E.A, 2); // Ap2 = A + 2
  fp2_sub_ui(&tmp[0], E.A, 2);
  fp2_neg(&tmp[0], tmp[0]);
  fp2_inv(&tmp[2], tmp[0]); // iAm2 = 1 / (2-A)
  fp2_add_ui(&tmp[0], E.A, 6);
  fp2_mul(&tmp[1], tmp[0], tmp[2]);
  fp2_mul(&(phiE->B), E.B, tmp[2]); // iB = B iAm2
  fp2_scalar_si(&(phiE->A), tmp[1], -2); // iA = -2 (A+6) iAm2
  MC_a24(&(phiE->A24), phiE->A);
}

/* Apply a 4-isogeny of Montgomery curves */
void iso4_eval_kummer(KP *Q,
		const iso4 iso,
		const KP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_mul(&tmp[0], P.x, P.z); // z1 = x z
  fp2_sub(&tmp[2], P.x, P.z);
  fp2_sqr(&tmp[1], tmp[2]); // x1 = (x - z)^2
  fp2_mul(&tmp[2], tmp[0], iso.Ap2); // zA2 = z1 Ap2
  fp2_scalar_si(&tmp[3], tmp[0], 4); // fourz = 4 z1
  fp2_add(&tmp[6], tmp[1], tmp[2]);
  fp2_add(&tmp[4], tmp[1], tmp[3]);
  fp2_mul(&tmp[5], tmp[6], tmp[4]); // x0 = (x1+zA2)(x1+fourz)

  fp2_sub(&tmp[4], tmp[3], tmp[2]);
  fp2_mul(&(Q->z), tmp[1], tmp[4]); // Z = x1 (fourz - zA2)
  fp2_copy(&(Q->x), tmp[5]);
}

/* Apply a 4-isogeny of Montgomery curves */
void iso4_eval_proj(MP *Q,
		const iso4 iso,
		const MP P) {
  fp2* tmp = params_global.fp2tmp;

  fp2_mul(&tmp[0], P.x, P.z); // z1 = x z
  fp2_sub(&tmp[2], P.x, P.z);
  fp2_sqr(&tmp[1], tmp[2]); // x1 = (x - z)^2
  fp2_mul(&tmp[2], tmp[0], iso.Ap2); // zA2 = z1 Ap2
  fp2_scalar_si(&tmp[3], tmp[0], 4); // fourz = 4 z1
  fp2_add(&tmp[6], tmp[1], tmp[2]);
  fp2_add(&tmp[4], tmp[1], tmp[3]);
  fp2_mul(&tmp[5], tmp[6], tmp[4]); // x0 = (x1+zA2)(x1+fourz)

  fp2_mul(&tmp[4], P.x, tmp[1]); // B = x x1
  fp2_mul(&tmp[5], tmp[5], tmp[4]); // x0 *= B
  fp2_sqr(&tmp[6], P.z);
  fp2_sub(&tmp[7], tmp[0], tmp[6]);
  fp2_scalar_si(&tmp[6], tmp[7], 2);
  fp2_add(&tmp[7], tmp[6], tmp[1]); // C = x1 + 2(z1 - z^2)
  fp2_sqr(&tmp[6], tmp[1]);
  fp2_mul(&tmp[8], tmp[2], tmp[3]);
  fp2_sub(&tmp[6], tmp[6], tmp[8]); // D = x1^2 - zA2 fourz
  fp2_mul(&tmp[8], tmp[7], tmp[6]);
  fp2_mul(&(Q->y), P.y, tmp[8]); // Y = y C D
  fp2_sqr(&tmp[6], tmp[4]);
  fp2_sub_ui(&tmp[7], iso.Ap2, 4);
      
  //The code below implements fp2_neg(&tmp[7], tmp[7]). For some reason just calling the function here has unpredictable behavious due to a bug in mpz_sub
  if (mpz_sgn(tmp[7].a) == 0)
    mpz_set(tmp[7].a, tmp[7].a);
  else
    mpz_sub(tmp[7].a, params_global.p, tmp[7].a); 
                                    
  if (mpz_sgn(tmp[7].b) == 0)
    mpz_set(tmp[7].b, P.x.b);
  else
    mpz_sub(tmp[7].b, params_global.p, tmp[7].b);
      
  fp2_mul(&tmp[8], tmp[6], tmp[7]);
  fp2_mul(&(Q->z), P.z, tmp[8]); // Z = z B^2 (4 - Ap2)
  fp2_copy(&(Q->x), tmp[5]);
}
 
/******* COMPOSITE ISOGENIES **************/

/* Implementation of a queue */

void PQ_copy(point_queue *res, const point_queue pq){
  KP_copy_point(&(res->point),pq.point);
  res->h = pq.h;
  res->next = pq.next;
  res->prev = pq.prev;
}

void PQ_print(point_queue pq, char * str){
  printf("\n%s\n",str);
  KP_print(pq.point,"point:");
  printf("\nh: %d", pq.h);
}

void PQ_init( point_queue *pq){
  KP_init(&(pq->point));
  pq->h=0;
  pq->next = malloc(sizeof(point_queue));
  pq->prev = malloc(sizeof(point_queue));
}

void PQ_clear(point_queue *pq){
  KP_clear(&(pq->point));
  free(pq->next);
  free(pq->prev);
}

/* Push (Px, Py, Pz) and (Qx, Qy, Qz) through the isogeny of kernel
 generated by (Rx, Rz) using the given strategy. */
void push_through_iso_alice(MC *E,
                      const MP R,
                      const int ell, int *strategy, int h,
                      MP *Pother,
                      MP *Qother, int e) {
	
  
  //fp2 *A = &(E->A); fp2 *B = &(E->B); fp2 *A24 = &(E->A24);
  int split, i, first = 1;
  union isogenies phi;
  point_queue *tail, *tmp;
  fp2_init(&phi.d1.u);
  fp2_init(&phi.d1.r);
  fp2_init(&phi.d2);
  fp2_init(&phi.d4.Ap2);
  Q_INIT(tail);
  KP_copy(&(tail->point),R.x,R.z);
  tail->h = h;
  while (tail) {
    h = tail->h;
    split = strategy[h];
    // Descend to the floor
    while (h > 1) {
      Q_INIT(tmp);
	  KP_copy(&(tmp->point),tail->point.x,tail->point.z);
      for ( i=0 ; i < h - split ; i++) {	  
      mont_double(&(tmp->point),
                   tmp->point, E->A24); 
      }
    tmp->h = split;
    Q_PUSH(tail, tmp);
    h = split;
    split = strategy[h];
    }
    // For ell=2, at the first iteration, bring the
    // 2-torsion point to (0,0)
    if (first) {
      first = 0;
      Q_INIT(tmp); // slight abuse
      mont_double(&(tmp->point), tail->point, E->A24);
      isom_comp(&phi.d1, E,
                *E, tmp->point);
      Q_CLEAR(tmp);
      APPLY_ISOG(isom_eval_kummer, isom_eval_proj, phi.d1, 0);
    }
    // Compute and apply the isogeny
    COMP_ISOG(iso2_comp, phi.d2);
    APPLY_ISOG(iso2_eval_kummer, iso2_eval_proj, phi.d2, 1);
  }
  // For ell=2 there is still a 4-isogeny to apply
  iso4_comp(&phi.d4, E, *E);
  // The queue is empty
  APPLY_ISOG(iso4_eval_kummer, iso4_eval_proj, phi.d4, 2);
  fp2_clear(&phi.d1.u);
  fp2_clear(&phi.d1.r);
  fp2_clear(&phi.d2);
  fp2_clear(&phi.d4.Ap2);
}

/* Push (Px, Py, Pz) and (Qx, Qy, Qz) through the isogeny of kernel
 generated by (Rx, Rz) using the given strategy. */
void push_through_iso_bob(MC *E,
                      const MP R,
                      const int ell, int *strategy, int h,
                      MP *Pother,
                      MP *Qother, int e) {
	
  
  //fp2 *A = &(E->A); fp2 *B = &(E->B); fp2 *A24 = &(E->A24);
  int split, i;
  union isogenies phi;
  point_queue *tail, *tmp;
  fp2_init(&phi.d3.p);
  fp2_init(&phi.d3.p2);
  Q_INIT(tail);
  KP_copy(&(tail->point),R.x,R.z);
  tail->h = h;
  while (tail) {
    h = tail->h;
    split = strategy[h];
    // Descend to the floor
    while (h > 1) {
      Q_INIT(tmp);
      KP_copy(&(tmp->point),tail->point.x,tail->point.z);
      for ( i=0 ; i < h - split ; i++) {	  
        mont_triple(&(tmp->point),
                     tmp->point, E->A24);
      }
      tmp->h = split;
    Q_PUSH(tail, tmp);
    h = split;
    split = strategy[h];
    }
    COMP_ISOG(iso3_comp, phi.d3);
    APPLY_ISOG(iso3_eval_kummer, iso3_eval_proj, phi.d3, 1);
  }
  fp2_clear(&phi.d3.p);
  fp2_clear(&phi.d3.p2);
}