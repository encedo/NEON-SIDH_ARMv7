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

#include "sidh.h"

void sidh_round1_alice(MP *Pother2, MP *Qother2, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, MP Pother, MP Qother, MP QP, MP PQ, int e){
  fp2 *A, *B, *A24;
  MP *Ptmp, *Qtmp, *PQtmp, *Rtmp;
  MC *E;

  Ptmp = malloc(sizeof(MP));
  MP_init(Ptmp);
  MP_copy(Ptmp,P);
  
  Qtmp = malloc(sizeof(MP));
  MP_init(Qtmp);
  MP_copy(Qtmp,Q);
  
  PQtmp = malloc(sizeof(MP));
  MP_init(PQtmp);
  MP_copy(PQtmp,PQ);
  
  Rtmp = malloc(sizeof(MP));
  MP_init(Rtmp);
  
  if(mpz_cmp_ui(M, 1)==0){
    MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
    MP_x_affine(PQtmp,QP);
	mont_3ladder(Rtmp, N, Ptmp->x, Qtmp->x, PQtmp->x, P.curve.A24); 
  }else if(mpz_cmp_ui(N, 1)==0){
	MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
	MP_x_affine(PQtmp,PQ);
	mont_3ladder(Rtmp, M, Qtmp->x, Ptmp->x, PQtmp->x, P.curve.A24);
  }else{
    shamir(Rtmp, *Ptmp, *Qtmp, M, N);  
  }

  E = malloc(sizeof(MC));
  MC_init(E);
  A = malloc(sizeof(fp2));
  fp2_init(A);
  fp2_copy(A,P.curve.A);
  B = malloc(sizeof(fp2));
  fp2_init(B);
  fp2_copy(B,P.curve.B);
  A24 = malloc(sizeof(fp2));
  fp2_init(A24);
  fp2_copy(A24,P.curve.A24);
    
  MC_set_curve(E, *A, *B, *A24);
  
  MP_copy_point(Ptmp,Pother);
  MP_copy_point(Qtmp,Qother);
  
  push_through_iso_alice( E, *Rtmp, ell, strat, len - 1, Ptmp, Qtmp, e);

  MP_set_curve(Pother2, *E);
  MP_set_curve(Qother2, *E);

  MP_copy_point(Pother2,*Ptmp);
  MP_copy_point(Qother2,*Qtmp);
    
  MP_clear(Ptmp); MP_clear(Qtmp); MP_clear(PQtmp); MP_clear(Rtmp);
  free(Ptmp); free(Qtmp); free(PQtmp); free(Rtmp);
  free(E);free(A);free(B); free(A24);
}

void sidh_round1_bob(MP *Pother2, MP *Qother2, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, MP Pother, MP Qother, MP QP, MP PQ, int e){
  fp2 *A, *B, *A24;
  MP *Ptmp, *Qtmp, *PQtmp, *Rtmp;
  MC *E;

  Ptmp = malloc(sizeof(MP));
  MP_init(Ptmp);
  MP_copy(Ptmp,P);
  
  Qtmp = malloc(sizeof(MP));
  MP_init(Qtmp);
  MP_copy(Qtmp,Q);
  
  PQtmp = malloc(sizeof(MP));
  MP_init(PQtmp);
  MP_copy(PQtmp,PQ);
  
  Rtmp = malloc(sizeof(MP));
  MP_init(Rtmp);
  
  if(mpz_cmp_ui(M, 1)==0){
    MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
    MP_x_affine(PQtmp,QP);
	mont_3ladder(Rtmp, N, Ptmp->x, Qtmp->x, PQtmp->x, P.curve.A24); 
  }else if(mpz_cmp_ui(N, 1)==0){
	MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
	MP_x_affine(PQtmp,PQ);
	mont_3ladder(Rtmp, M, Qtmp->x, Ptmp->x, PQtmp->x, P.curve.A24);
  }else{
    shamir(Rtmp, *Ptmp, *Qtmp, M, N);  
  }

  E = malloc(sizeof(MC));
  MC_init(E);
  A = malloc(sizeof(fp2));
  fp2_init(A);
  fp2_copy(A,P.curve.A);
  B = malloc(sizeof(fp2));
  fp2_init(B);
  fp2_copy(B,P.curve.B);
  A24 = malloc(sizeof(fp2));
  fp2_init(A24);
  fp2_copy(A24,P.curve.A24);
    
  MC_set_curve(E, *A, *B, *A24);
  
  MP_copy_point(Ptmp,Pother);
  MP_copy_point(Qtmp,Qother);
  
  push_through_iso_bob( E, *Rtmp, ell, strat, len - 1, Ptmp, Qtmp, e);

  MP_set_curve(Pother2, *E);
  MP_set_curve(Qother2, *E);

  MP_copy_point(Pother2,*Ptmp);
  MP_copy_point(Qother2,*Qtmp);
    
  MP_clear(Ptmp); MP_clear(Qtmp); MP_clear(PQtmp); MP_clear(Rtmp);
  free(Ptmp); free(Qtmp); free(PQtmp); free(Rtmp);
  free(E);free(A);free(B); free(A24);
}

void sidh_round2_alice(MC *E, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, int e){ 
	
  fp2 *A, *B, *A24;
  MP *Ptmp, *Qtmp, *PQtmp, *Rtmp;

  Ptmp = malloc(sizeof(MP));
  MP_init(Ptmp);
  MP_copy(Ptmp,P);
  
  Qtmp = malloc(sizeof(MP));
  MP_init(Qtmp);
  MP_copy(Qtmp,Q);
  
  PQtmp = malloc(sizeof(MP));
  MP_init(PQtmp);
  MP_copy(PQtmp,P);
  
  Rtmp = malloc(sizeof(MP));
  MP_init(Rtmp);
    
  if(mpz_cmp_ui(M, 1)==0){
	MP_sub(PQtmp,Q,P);
    MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
    MP_x_affine(PQtmp,*PQtmp);
	mont_3ladder(Rtmp, N, Ptmp->x, Qtmp->x, PQtmp->x, P.curve.A24); 
  }else if(mpz_cmp_ui(N, 1)==0){
	MP_sub(PQtmp,P,Q);
	MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
	MP_x_affine(PQtmp,*PQtmp);
	mont_3ladder(Rtmp, M, Qtmp->x, Ptmp->x, PQtmp->x, P.curve.A24);
  }else{  
    shamir(Rtmp, *Ptmp, *Qtmp, M, N);  
  }
  
  A = malloc(sizeof(fp2));
  fp2_init(A);
  fp2_copy(A,P.curve.A);
  B = malloc(sizeof(fp2));
  fp2_init(B);
  fp2_copy(B,P.curve.B);  
  A24 = malloc(sizeof(fp2));
  fp2_init(A24);
  fp2_copy(A24,P.curve.A24);
  
  MC_set_curve(E, *A, *B, *A24);
  push_through_iso_alice( E, *Rtmp, ell, strat, len - 1, NULL, NULL, e);
  
  
  MP_clear(Ptmp); MP_clear(Qtmp); MP_clear(PQtmp); MP_clear(Rtmp);
  free(Ptmp); free(Qtmp); free(PQtmp); free(Rtmp);
  free(A);free(B); free(A24);
}

void sidh_round2_bob(MC *E, mpz_t M, mpz_t N, int ell, int *strat, int len, MP P, MP Q, int e){ 
	
  fp2 *A, *B, *A24;
  MP *Ptmp, *Qtmp, *PQtmp, *Rtmp;

  Ptmp = malloc(sizeof(MP));
  MP_init(Ptmp);
  MP_copy(Ptmp,P);
  
  Qtmp = malloc(sizeof(MP));
  MP_init(Qtmp);
  MP_copy(Qtmp,Q);
  
  PQtmp = malloc(sizeof(MP));
  MP_init(PQtmp);
  MP_copy(PQtmp,P);
  
  Rtmp = malloc(sizeof(MP));
  MP_init(Rtmp);
    
  if(mpz_cmp_ui(M, 1)==0){
	MP_sub(PQtmp,Q,P);
    MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
    MP_x_affine(PQtmp,*PQtmp);
	mont_3ladder(Rtmp, N, Ptmp->x, Qtmp->x, PQtmp->x, P.curve.A24); 
  }else if(mpz_cmp_ui(N, 1)==0){
	MP_sub(PQtmp,P,Q);
	MP_x_affine(Ptmp,P);
	MP_x_affine(Qtmp,Q);
	MP_x_affine(PQtmp,*PQtmp);
	mont_3ladder(Rtmp, M, Qtmp->x, Ptmp->x, PQtmp->x, P.curve.A24);
  }else{  
    shamir(Rtmp, *Ptmp, *Qtmp, M, N);  
  }
  
  A = malloc(sizeof(fp2));
  fp2_init(A);
  fp2_copy(A,P.curve.A);
  B = malloc(sizeof(fp2));
  fp2_init(B);
  fp2_copy(B,P.curve.B);  
  A24 = malloc(sizeof(fp2));
  fp2_init(A24);
  fp2_copy(A24,P.curve.A24);
  
  MC_set_curve(E, *A, *B, *A24);
  push_through_iso_bob( E, *Rtmp, ell, strat, len - 1, NULL, NULL, e);
  
  
  MP_clear(Ptmp); MP_clear(Qtmp); MP_clear(PQtmp); MP_clear(Rtmp);
  free(Ptmp); free(Qtmp); free(PQtmp); free(Rtmp);
  free(A);free(B); free(A24);
}


   // returns 0 if shared keys are different, 1 if they are the same
int sidh_full(double *timeKey, double *timeAliceA, double *timeBobA, double *timeAliceB, double* timeBobB, char * eA, char * eB, char * lA_str, char * lB_str, int *strA, int lenA, int *strB, int lenB, MP *PA, MP *QA, MP *PB, MP *QB){
  mpz_t *mA, *nA, *mB, *nB;
  struct timeval tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8, diff1, diff2, diff3, diff4, total;
  int lA, lB;
  int good=0;
  
  mA = malloc(sizeof(mpz_t));
  nA = malloc(sizeof(mpz_t));
  mB = malloc(sizeof(mpz_t));
  nB = malloc(sizeof(mpz_t));
  mpz_init(mA);
  mpz_init(nA);
  mpz_init(mB);
  mpz_init(nB);
   
  if (VERBOSE) printf("\n\n----------------Randomly Generating Secret Keys\n");
  rand_subgroup(mA,nA,lA_str,eA);
  rand_subgroup(mB,nB,lB_str,eB);
  
  /*
  mpz_set_str(mA,"33870763417823928285299756663564920023337704159476068098191993835048334605528",10);
  mpz_set_str(nA,"1",10);
  mpz_set_str(mB,"162610922938997627305219831441371930844092661195644947061471713683209124997675",10);
  mpz_set_str(nB,"1",10);
  */
  
  lA = atoi(lA_str);
  lB = atoi(lB_str);
  
  MC EsA, EsB;
  MP phiPB, phiQB, phiPA, phiQA, QPB, PQB, QPA, PQA ;
  MP_init(&phiPB);
  MP_init(&phiQB);
  MP_init(&phiPA);
  MP_init(&phiQA);
    
  MC_init(&EsA);
  MC_init(&EsB);
    
  MP_init(&QPB);
  MP_init(&PQB);
  MP_init(&QPA);
  MP_init(&PQA);
    
  MP_sub(&QPB, *QB, *PB);
  MP_sub(&PQB, *PB, *QB);
  MP_sub(&QPA, *QA, *PA);
  MP_sub(&PQA, *PA, *QA);
    
  int eA_num = atoi(eA);
  int eB_num = atoi(eB);
  if (VERBOSE) printf("----------------Alice1 Generating Alice's Ephemeral Key\n");
  gettimeofday(&tv1,NULL);
  sidh_round1_alice(&phiPB, &phiQB, *mA, *nA, lA, strA, lenA, *PA, *QA, *PB, *QB, QPA, PQA, eA_num);
  gettimeofday(&tv2,NULL);
  if (VERBOSE) printf("----------------Bob1 Generating Bob's Ephemeral Key\n");
  gettimeofday(&tv3,NULL);
  sidh_round1_bob(&phiPA, &phiQA, *mB, *nB, lB, strB, lenB, *PB, *QB, *PA, *QA, QPB, PQB, eB_num);
  gettimeofday(&tv4,NULL);
  if (VERBOSE) printf("----------------Alice2 Computing Shared Key on Alice's Side\n");
  gettimeofday(&tv5,NULL);
  sidh_round2_alice(&EsA, *mA, *nA, lA, strA, lenA, phiPA, phiQA, eA_num);
  gettimeofday(&tv6,NULL);
  if (VERBOSE) printf("----------------Bob2 Computing Shared Key on Bob's Side\n");
  gettimeofday(&tv7,NULL);
  sidh_round2_bob(&EsB, *mB, *nB, lB, strB, lenB, phiPB, phiQB, eB_num);
  gettimeofday(&tv8,NULL);
    
  fp2 EsA_j, EsB_j;
  fp2_init( &EsA_j);
  fp2_init( &EsB_j);
    
  MC_j_invariant(&EsA_j, EsA);
  MC_j_invariant(&EsB_j, EsB);
  if (VERBOSE) printf("\n\n");
    
  if ( fp2_equals(&EsA_j, &EsB_j)!=1 ){
	if (VERBOSE){
      gmp_printf("ERROR: the shared keys don't match! Here's the secret keys:\n\tmA = %Zd\n\tnA = %Zd\n\tmB = %Zd\n\tnB = %Zd\n",*mA,*nA,*mB,*nB);
        
      printf("*************EsA*************\n");
      MC_print(EsA);  
      fp2_print(EsA_j,"EsA J-Invariant");
      printf("\n*************EsB*************\n");
      MC_print(EsB); 
      fp2_print(EsB_j,"EsB J-Invariant");	  
	}
    good = 0;
  }else{
	if (VERBOSE){
      printf("*************Shared Key(s)*************\n\n");
      printf("*************EsA*************\n");
      MC_print(EsA);  
      fp2_print(EsA_j,"EsA J-Invariant");
      printf("\n*************EsB*************\n");
      MC_print(EsB); 
      fp2_print(EsB_j,"EsB J-Invariant");
      printf("\n*************Secret Keys*************\n");
      gmp_printf(" mA: %Zd\n nA: %Zd\n mB: %Zd\n nB: %Zd\n",*mA,*nA,*mB,*nB);    
      printf("\n*************J-Invariants Match*************\n\n");	  
	}
	good = 1;
  }
  free(mA);free(mB);free(nA);free(nB);
  MP_clear(&phiPB); MP_clear(&phiQB); MP_clear(&phiPA); MP_clear(&phiQA);
  MC_clear(&EsA); MC_clear(&EsB);
  MP_clear(&QPB); MP_clear(&PQB); MP_clear(&QPA); MP_clear(&PQA);
  fp2_clear(&EsA_j); fp2_clear(&EsB_j);
   
  timersub(&tv2, &tv1, &diff1);
  timersub(&tv4, &tv3, &diff2);
  timersub(&tv6, &tv5, &diff3);
  timersub(&tv8, &tv7, &diff4);
    
  timeradd(&diff1, &diff2, &total);
  timeradd(&diff3, &total, &total);
  timeradd(&diff4, &total, &total);
  
  if (VERBOSE){
    printf ("Keygen 1 Real Time (sec) :%f \n", (double) diff1.tv_usec/1000000 + diff1.tv_sec);
    printf ("Keygen 2 Real Time (sec) :%f \n", (double) diff2.tv_usec/1000000 + diff2.tv_sec);
    printf ("Keygen 3 Real Time (sec) :%f \n", (double) diff3.tv_usec/1000000 + diff3.tv_sec);
    printf ("Keygen 4 Real Time (sec) :%f \n", (double) diff4.tv_usec/1000000 + diff3.tv_sec);
    printf ("Total Time (sec) :%f \n", (double) total.tv_usec/1000000.0 + total.tv_sec);	
  }

    
  *timeKey = (double) total.tv_usec/1000000.0 + total.tv_sec;
  *timeAliceA = (double) diff1.tv_usec/1000000.0 + diff1.tv_sec;
  *timeBobA = (double) diff2.tv_usec/1000000.0 + diff2.tv_sec;
  *timeAliceB = (double) diff3.tv_usec/1000000.0 + diff3.tv_sec;
  *timeBobB = (double) diff4.tv_usec/1000000.0 + diff4.tv_sec;
  
  return good;
}