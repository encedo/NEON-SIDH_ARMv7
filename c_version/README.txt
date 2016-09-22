NEON-SIDH v1.0 (C version)

The NEON-SIDH library is a supersingular isogeny-based cryptography library that implements
the supersingular isogeny Diffie-Hellman (SIDH) key exchange protocol. The current version of the library
targets the affine isogeny formulas provided by Luca de Feo, David Jao, and Jerome Plut's paper, "Towards Quantum-resistant Cryptosystems 
from Supersingular Elliptic Curve Isogenies".

This library was developed to test how efficient the SIDH protocol can be on ARMv7 devices. This is detailed
in the accompanying paper:

NEON-SIDH: Efficient Implementation of Supersingular Isogeny Diffie-Hellman
Key Exchange Protocol on ARM

Authors: Brian Koziel, Amir Jalali, Reza Azarderakhsh, David Jao, and
	Mehran Mozaffari-Kermani
	
Link to published paper: To be posted 


This is the C version of this protocol that is based on the GNU Multiprecision Library (GMP) for standard
arithmetic. GMP features extremely fast, but not constant-time, arithmetic. Notably, the use of this library
allows for portable optimizations between platforms. This code is intended to be developed in a Linux environment
(GCC) running on an ARMv7 CPU and has been verified on GMP v6.1.0.

To compile on Linux using GNU GCC, execute the following command from the command prompt:

make ARM_CORE=CORTEX-A[5,7,8,9,12,15,17]

This generates an executable, sidh_test, that simulates and times a full key exchange using SIDH. 
This can be run with the following command:

./sidh_test <parameter_file> <number of iterations>

The parameter file must conform to a very standard form, shown in parameter_file_order.txt. 
Several parameter files for various prime sizes and curves are provided, but more can be 
generated from Luca de Feo's Sage code at https://github.com/defeo/ss-isogeny-software . This code 
attempts to work with any valid sets of parameters and private keys.

