# OMP-DGHV
This is a toy example of a homomorphic encryption scheme, namely the Fully Homomorphic Encryption scheme over the Integers augmented with flattening techniques and enhanced with omp support.

This code is an implementation of the DGHV scheme as described in "Fully Homomorphic Encryption over the Integers"
by van Dijk et.al., available at https://eprint.iacr.org/2009/616.pdf, altogheter with flattening techniques described by 
Gentry et.al. in "Homomorphic Encryption from Learning with Errors : Conceptually-Simpler, Asymptotically-Faster, Attribute-Based",
https://eprint.iacr.org/2013/340.pdf. More exaxtly, the implementation is a result of one the suggestions found in the apendix
of the later paper.

This work was also inspired by "Flattening NTRU for Evaluation Key Free Homomorphic Encryption", by Doroz and Sunar, which can
be found at https://eprint.iacr.org/2016/315.pdf.

System requirements :
1. [NTL](http://www.shoup.net/ntl/): A Library for doing Number Theory 9.5.0 (requires C++11) NOTE: to avoid random crashes compile it running
./configure NTL_EXCEPTIONS=o linked with

2. [gmp library version gmp-6.1.0.](https://gmplib.org/)

3. [openmp 3.0](http://openmp.org/wp/), if speedup is desired and the hardware system is able to support multithreading
The code works without openmp, but it is slower :).

4. a compiler with support for [c++11](https://en.wikipedia.org/wiki/C%2B%2B11)
