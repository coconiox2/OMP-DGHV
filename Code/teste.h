// DGHV securiy parameter
// dupa Fully homomorphic encryption over the Integers with Shorter Public Key by Coron et.al
// http://publications.uni.lu/bitstream/10993/12396/1/441.pdf
// #define SEC_PARAM 42		// toy 
// #define SEC_PARAM 52		// small
// #define SEC_PARAM 62		// medium 
// #define SEC_PARAM 72		// large

# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <omp.h>
# include <time.h>

#include "utilities.h"
#include <time.h>
#include <assert.h>
#include "Flat_DGHV.h"
#include <fstream>
#include <NTL/ZZ.h>

#include "../Applications/HE_Utils.h"
#include "../Applications/HE_Integer.h"

#define SEC_PARAM 140
#define NR_TESTE 1000
// #define MULT_DEPTH 50

using namespace std;

extern Flat_DGHV* HE_Context;
extern Mat_ZZ* ctxt_of_1;                  // encryption of 1 (constant)
extern int		t_bits;
extern ZZ *baza;

void test_max_func(int lambda = 192, int no_bits = 8, int N = 16);

void test_FDGHV_max_depth();

void test_Flat_DGHV_simple();

vector<int> compute_random_circuit(int size = NR_TESTE );

void test_schema_FDGHV();

void test_max_mult_depth();

void test_max_no_encryptions();

void cmp_HElib();








