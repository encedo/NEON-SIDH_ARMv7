####  Makefile for compilation on Linux  ####

#	NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
#		Key Exchange Protocol on ARM
#		
#	Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
#		Mehran Mozaffari-Kermani
#		
#	Link to published paper: To be posted 
#		
#	Previous versions by Luca de Feo and Dieter Fishbein
#	
#	ecc.h: Elliptic curve cryptography (ECC) functions, primarily
#		focused on Montgomery curves
#	
#	Copyright (c) Brian Koziel 2016
#
#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


OPT=-O3 -fomit-frame-pointer -flto -fprefetch-loop-arrays -funroll-all-loops -ftree-vectorize -ftree-vectorizer-verbose=0

COMPILER=gcc

ARM_FPU=-mfpu=neon

ifeq "$(ARM_CORE)" "CORTEX-A5"
    ARCHITECTURE=-mtune=cortex-a5
else ifeq "$(ARM_CORE)" "CORTEX-A7"
    ARCHITECTURE=-mtune=cortex-a7
else ifeq "$(ARM_CORE)" "CORTEX-A8"
    ARCHITECTURE=-mtune=cortex-a8
else ifeq "$(ARM_CORE)" "CORTEX-A9"
    ARCHITECTURE=-mtune=cortex-a9
else ifeq "$(ARM_CORE)" "CORTEX-A12"
    ARCHITECTURE=-mtune=cortex-a12
else ifeq "$(ARM_CORE)" "CORTEX-A15"
    ARCHITECTURE=-mtune=cortex-a15
else ifeq "$(ARM_CORE)" "CORTEX-A17"
    ARCHITECTURE=-mtune=cortex-a17
endif

CC=$(COMPILER)
CFLAGS=-c
ARMFLAGS= $(OPT) $(ARM_FPU) $(ARCHITECTURE)
LIBFLAGS=-lm -w -lpthread -lrt -lgmp -lcrypto

OBJECTS=fp.o fp2.o ecc.o mont_curve.o iso.o sidh.o util.o
OBJECTS_SIDH_TEST=sidh_test.o
OBJECTS_ALL=$(OBJECTS) $(OBJECTS_SIDH_TEST)

sidh_test: $(OBJECTS_ALL)
	$(CC) $(ARMFLAGS) $(OBJECTS_ALL) -o sidh_test $(LIBFLAGS)

fp.o: fp.c fp.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) fp.c $(LIBFLAGS)
	
fp2.o: fp2.c fp2.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) fp2.c $(LIBFLAGS)
	
ecc.o: ecc.c ecc.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) ecc.c $(LIBFLAGS)
	
mont_curve.o: mont_curve.c mont_curve.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) mont_curve.c $(LIBFLAGS)
	
iso.o: iso.c iso.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) iso.c $(LIBFLAGS)
	
sidh.o: sidh.c sidh.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) sidh.c $(LIBFLAGS)
	
util.o: util.c util.h sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) util.c $(LIBFLAGS)

sidh_test.o: sidh_test.c sidh_params.h
	$(CC) $(CFLAGS) $(ARMFLAGS) sidh_test.c $(LIBFLAGS)

clean:
	rm sidh_test $(OBJECTS_ALL)

