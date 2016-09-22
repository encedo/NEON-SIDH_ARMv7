/*
	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
		Key Exchange Protocol on ARM
		
	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
		Mehran Mozaffari-Kermani
		
	Link to published paper: To be posted 
		
	Previous versions by Luca de Feo and Dieter Fishbein
	
	sidh_test.c: Full setup and simulation of SIDH protocol
	
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

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <openssl/rand.h>
#include <sys/time.h>
#include "ecc.h"
#include "fp.h"
#include "fp2.h"
#include "mont_curve.h"
#include "util.h"
#define DEBUG
#define CLKFREQ 2300

//1st argument specifies file with parameters, second # of times to run the key exchange. If second argument is null, only 1 iteration is performed
int main(int argc, char *argv[]) {
  int iterations;
  int i;
    
  if( argc < 2 ){
    printf("ERROR: Must specify file name in 1st command line agrument.\n");    
  }else{
    printf("File where Public parameters are: '%s'\n", argv[1]);
  }
    
  if( !argv[2] )
    iterations = 1;
  else{
    iterations = atoi(argv[2]);
    printf("Number of Iterations: %s\n", argv[2]);
  }    
    
  int MAX_LENGTH = 10000;
  int *strA, *strA_t, *strB, *strB_t; 
  int lenA=0;
  int lenB=0;
  MP *PA;
  MP *QA;
  MP *PB;
  MP *QB;
  char *p, *eA, *eB, *lA, *lB, *fMult;
    
  p = malloc(sizeof(char)*MAX_LENGTH);
  eA = malloc(sizeof(char)*MAX_LENGTH);
  eB = malloc(sizeof(char)*MAX_LENGTH);
  lA = malloc(sizeof(char)*MAX_LENGTH);
  lB = malloc(sizeof(char)*MAX_LENGTH);
  fMult = malloc(sizeof(char)*MAX_LENGTH);
  PA = malloc(sizeof(MP));
  QA = malloc(sizeof(MP));
  PB = malloc(sizeof(MP));
  QB = malloc(sizeof(MP));
  strA_t = malloc(MAX_LENGTH*sizeof(int));
  strB_t = malloc(MAX_LENGTH*sizeof(int));
    
  params_from_file( p, eA, eB, lA, lB, fMult, strA_t, &lenA, strB_t, &lenB, PA, QA, PB, QB, argv[1] );

  strA = malloc(lenA*sizeof(int));
  strB = malloc(lenB*sizeof(int));
    
  //Putting strategies into arrays of the exact required length
  int r=0;
  for(r; r<lenA; r++)
    strA[r] = strA_t[r];    
  int k=0;
  for(k; k<lenB; k++)
    strB[k] = strB_t[k]; 
	
  //run the key exchange 'iterations' number of times and calculate the arverage time taken
  i=0;
  int good = 0;
  double *timesKey = malloc(iterations*sizeof(double));
  double *timesAliceA = malloc(iterations*sizeof(double));
  double *timesBobA = malloc(iterations*sizeof(double));
  double *timesAliceB = malloc(iterations*sizeof(double));
  double *timesBobB = malloc(iterations*sizeof(double));
  double totalKey=0;
  double totalAliceA = 0;
  double totalBobA = 0;
  double totalAliceB = 0;
  double totalBobB = 0;
  double avgKey, medianKey, maxKey, minKey;
  double avgAliceA, medianAliceA, maxAliceA, minAliceA;
  double avgBobA, medianBobA, maxBobA, minBobA;
  double avgAliceB, medianAliceB, maxAliceB, minAliceB;
  double avgBobB, medianBobB, maxBobB, minBobB;
  int errors=0;
	
  maxKey = 0.0; minKey = 100.0;
  maxAliceA = 0.0; maxBobA = 0.0; maxAliceB = 0.0; maxBobB = 0.0;
  minAliceA = 100.0; minBobA = 100.0; minAliceB = 100.0; minBobB = 100.0; 
	
  double *timeKey;
  double *timeAliceA, *timeBobA, *timeAliceB, *timeBobB;
	
  timeKey = malloc(sizeof(double));
  timeAliceA = malloc(sizeof(double));
  timeBobA = malloc(sizeof(double));
  timeAliceB = malloc(sizeof(double));
  timeBobB = malloc(sizeof(double));
    
  for(i; i<iterations; i++){ 
    good = sidh_full(timeKey,timeAliceA,timeBobA,timeAliceB,timeBobB,eA,eB,lA,lB,strA,lenA,strB,lenB,PA,QA,PB,QB);
    totalKey += *timeKey;
    totalAliceA += *timeAliceA;
    totalBobA += *timeBobA;
    totalAliceB += *timeAliceB;
    totalBobB += *timeBobB;
	
    if (*timeKey > maxKey) maxKey = *timeKey;
    if (*timeAliceA > maxAliceA) maxAliceA = *timeAliceA;
    if (*timeBobA > maxBobA) maxBobA = *timeBobA;
	if (*timeAliceB > maxAliceB) maxAliceB = *timeAliceB;
    if (*timeBobB > maxBobB) maxBobB = *timeBobB;
	
    if (*timeKey < minKey) minKey = *timeKey;
    if (*timeAliceA < minAliceA) minAliceA = *timeAliceA;
    if (*timeBobA < minBobA) minBobA = *timeBobA;
    if (*timeAliceB < minAliceB) minAliceB = *timeAliceB;
    if (*timeBobB < minBobB) minBobB = *timeBobB;
	
    timesKey[i] = *timeKey;
    timesAliceA[i] = *timeAliceA;
    timesBobA[i] = *timeBobA;
    timesAliceB[i] = *timeAliceB;
    timesBobB[i] = *timeBobB;
	
    //printf("\nDecomp Total time: %f",totalTime2);
    printf ("\n                            Cur            Ave            Med             Max             Min\n");
    printf ("Key Exchange Time    (sec) :%f       %f       %f        %f        %f \n", *timeKey,totalKey/((i+1)*1.0),getMedian(timesKey,i+1),maxKey,minKey);
    printf ("Alice (A)    Time    (sec) :%f       %f       %f        %f        %f \n", *timeAliceA,totalAliceA/((i+1)*1.0),getMedian(timesAliceA,i+1),maxAliceA,minAliceA);
    printf ("Bob   (A)    Time    (sec) :%f       %f       %f        %f        %f \n", *timeBobA,totalBobA/((i+1)*1.0),getMedian(timesBobA,i+1),maxBobA,minBobA);
    printf ("Alice (B)    Time    (sec) :%f       %f       %f        %f        %f \n", *timeAliceB,totalAliceB/((i+1)*1.0),getMedian(timesAliceB,i+1),maxAliceB,minAliceB);
    printf ("Bob   (B)    Time    (sec) :%f       %f       %f        %f        %f \n", *timeBobB,totalBobB/((i+1)*1.0),getMedian(timesBobB,i+1),maxBobB,minBobB);
	
    if (!good)
      errors +=1;
    printf("\nErrors/Total so far: %d/%d\n", errors,i+1);
  }
    
  avgKey = totalKey/(iterations*1.0);
  avgAliceA = totalAliceA/(iterations*1.0);
  avgBobA = totalBobA/(iterations*1.0);
  avgAliceB = totalAliceB/(iterations*1.0);
  avgBobB = totalBobB/(iterations*1.0);
	
  medianKey = getMedian(timesKey,iterations);
  medianAliceA = getMedian(timesAliceA,iterations);
  medianBobA = getMedian(timesBobA,iterations);
  medianAliceB = getMedian(timesAliceB,iterations);
  medianBobB = getMedian(timesBobB,iterations);
	
  printf("\nApplication, Average,   Median,   Median Cycles for %d iteration(s) (sec)", iterations);
  printf("\nKey    :     %f\t%f\t%f",avgKey,medianKey, medianKey*CLKFREQ);
  printf("\nAliceA :     %f\t%f\t%f",avgAliceA,medianAliceA, medianAliceA*CLKFREQ);
  printf("\nBobA   :     %f\t%f\t%f",avgBobA,medianBobA, medianBobA*CLKFREQ);
  printf("\nAliceB :     %f\t%f\t%f",avgAliceB,medianAliceB, medianAliceB*CLKFREQ);
  printf("\nBobB   :     %f\t%f\t%f",avgBobB,medianBobB, medianBobB*CLKFREQ);
  printf("\nAliceT :     %f\t%f\t%f",avgAliceA+avgAliceB,medianAliceA+medianAliceB, (medianAliceA+medianAliceB)*CLKFREQ);
  printf("\nBobT   :     %f\t%f\t%f",avgBobA+avgBobB,medianBobA+medianBobB, (medianBobA+medianBobB)*CLKFREQ);
  
  printf("\nParameters: %s\n",  argv[1]);

  MC_clear(&PA->curve);
  MP_clear(PA); MP_clear(QA); MP_clear(PB); MP_clear(QB);
  
  free(timesKey); free(timesAliceA); free(timesBobA); free(timesAliceB); free(timesBobB);
	
  free(timeKey); free(timeAliceA); free(timeBobA); free(timeAliceB); free(timeBobB);
	
  free(p); free(eA); free(eB); free(lA); free(lB); free(fMult); free(PA); free(QA); 
  free(PB); free(QB); free(strA_t); free(strB_t); free(strA); free(strB);
}