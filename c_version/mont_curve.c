/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	mont_curve.c: Montgomery curve and point declarations and methods
	
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

#include "mont_curve.h"

/* Utility routine to compute (A+2)/4 */
void MC_a24(fp2* A24, const fp2 A) {
  fp2* tmp = params_global.fp2tmp;

  fp2_add_ui(&tmp[0], A, 2);
  fp2_scalar(A24, tmp[0], params_global.four_inverse);
}


void MC_copy(MC *res, const MC curve){
  fp2_copy( &(*res).A, curve.A );
  fp2_copy( &(*res).B, curve.B );
  fp2_copy( &(*res).A24, curve.A24 );
}

void MC_init(MC* x) {
  fp2_init( &(*x).A);
  fp2_init( &(*x).B);
  fp2_init( &(*x).A24);
}

void MC_clear(MC *x) {
  fp2_clear(&(*x).A);
  fp2_clear(&(*x).B);
  fp2_clear(&(*x).A24);
}

void MC_print(const MC x){
  fp2_print( x.A,"A" );
  fp2_print( x.B,"B" );
  fp2_print( x.A24,"A24" ); 
}

void MC_set_curve( MC *curve, fp2 A, fp2 B, fp2 A24 ){
  (*curve).A = A;
  (*curve).B = B;
  (*curve).A24 = A24;
}

void MC_j_invariant(fp2 *final, const MC curve){
  fp2 tmp1, tmp2, tmp3, denom, num;
    
  fp2_init(&tmp1);
  fp2_init(&tmp2);
  fp2_init(&tmp3);
  fp2_init(&denom);
  fp2_init(&num);
  fp2_sqr(&tmp1, curve.A);
  fp2_sub_ui(&denom, tmp1, 4);
  fp2_sub_ui(&tmp2, tmp1, 3 );
  fp2_sqr( &tmp3, tmp2);
  fp2_mul( &tmp1, tmp3, tmp2);
  fp2_scalar_si(&num, tmp1, 256);
  fp2_div(final, num, denom);
  
  fp2_clear(&tmp1); fp2_clear(&tmp2); fp2_clear(&tmp3); fp2_clear(&denom); fp2_clear(&num);
}

void MP_init(MP* a) {
  fp2_init( &(a->x));
  fp2_init( &(a->y));
  fp2_init( &(a->z));
  MC_init(&(a->curve));
}

void MP_copy_point(MP *res, const MP P){
  fp2_copy(&(res->x),P.x);
  fp2_copy(&(res->y),P.y);
  fp2_copy(&(res->z),P.z);  
}
    
void MP_copy( MP *res, const MP P ){
  fp2_copy(&(res->x),P.x);
  fp2_copy(&(res->y),P.y);
  fp2_copy(&(res->z),P.z);  
  MC_copy(&(res->curve),P.curve);
}
    
void MP_set_curve( MP *res, const MC curve ){
  (*res).curve.A = curve.A;
  (*res).curve.B = curve.B;
  (*res).curve.A24 = curve.A24;
}

void MP_set( MP *res, const fp2 x, const fp2 y, const fp2 z, const MC curve ){
  (*res).x = x;
  (*res).y = y;
  (*res).z = z;
    
  (*res).curve.A = curve.A;
  (*res).curve.B = curve.B;
  (*res).curve.A24 = curve.A24;
}

void MP_clear(MP *a) {
  fp2_clear(&(*a).x);
  fp2_clear(&(*a).y);
  fp2_clear(&(*a).z);
}

void MP_print(const MP a, char * string){
  printf("%s: \n",string);
  fp2_print( a.x,"x" );
  fp2_print( a.y,"y" );
  fp2_print( a.z,"z" );
  MC_print(  a.curve );
}

//P and Q must be on the same curve
void MP_add(MP *res, const MP P, const MP Q){
  fp2* tmp = params_global.fp2tmp;
    
  fp2_mul( &tmp[1], P.y, Q.z );
  fp2_mul( &tmp[2] , P.x, Q.z );
  fp2_mul( &tmp[3] , P.z, Q.z );
  fp2_mul( &tmp[4], Q.y ,P.z );
    
  fp2_sub( &tmp[5], tmp[4], tmp[1]);
  fp2_sqr( &tmp[6], tmp[5] );
  fp2_mul( &tmp[4], Q.x, P.z );
  fp2_sub( &tmp[7], tmp[4], tmp[2] );
	
  fp2_sqr( &tmp[8], tmp[7] );
  fp2_mul( &tmp[0], tmp[7], tmp[8] );
  fp2_mul( &tmp[2], tmp[8], tmp[2] );

  fp2_mul( &tmp[4], P.curve.B, tmp[6] );
  fp2_mul( &tmp[1], P.curve.A, tmp[8] );
  fp2_sub( &tmp[8], tmp[4], tmp[1] );
  fp2_mul( &tmp[1], tmp[8], tmp[3] );
  fp2_sub( &tmp[4], tmp[1], tmp[0] );
  fp2_scalar_si( &tmp[1], tmp[2], 2 );
  fp2_sub( &tmp[4], tmp[4], tmp[1] ); 

  fp2_mul(&(*res).x, tmp[7], tmp[4] );
  fp2_mul(&(*res).z, tmp[0], tmp[3] );

  fp2_sub( &tmp[4], tmp[2], tmp[4] );
  fp2_mul( &tmp[1], tmp[5], tmp[4] );   
  fp2_mul( &tmp[6], tmp[0], tmp[1] );
  fp2_sub(&(*res).y, tmp[1], tmp[6] );
  (*res).curve = P.curve;
}

//P and Q must be on the same curve
void  MP_sub(MP *res, const MP P, const MP Q){
  MP *point;
  point = malloc(sizeof(MP));
  MP_init(point);
    
  MP_neg(point,Q);
  MP_add(res,P,*point);
  free(point);
}

//Returns the negative of P
void MP_neg(MP *res, const MP P){
  fp2 *yN;
    
  yN = malloc(sizeof(fp2));
  fp2_init( yN);
  fp2_neg( yN , (P).y );

  (*res).x = (P).x;
  (*res).y = *yN;
  (*res).z = (P).z;
  (*res).curve = P.curve;
  free(yN);
}

void MP_x_affine(MP *res, const MP P){
  fp2_div(&(res->x), (res->x),(res->z));
  fp2_set(&(res->z),"0","1");
  //(*res).y = *yN; //Don't change y
  //(*res).z = (P).z;
  res->curve = P.curve;
}

void KP_init(KP* a){
  fp2_init(&(a->x));
  fp2_init(&(a->z));
}

void KP_clear(KP *a){
  fp2_clear(&(a->x));
  fp2_clear(&(a->z));
}

void KP_copy(KP *res, const fp2 x, const fp2 z){
  fp2_copy(&(res->x),x);
  fp2_copy(&(res->z),z);
}
    
void KP_copy_point(KP *res, const KP P){
  fp2_copy(&(res->x),P.x);
  fp2_copy(&(res->z),P.z);   
}

void KP_set( KP *res, const fp2 x, const fp2 z){
  (*res).x = x;
  (*res).z = z;
}

void KP_print(const KP a, char * string){
  printf("%s: \n",string);
  fp2_print( a.x,"x" );
  fp2_print( a.z,"z" );  
}