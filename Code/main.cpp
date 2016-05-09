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

// DGHV securiy parameter
// dupa Fully homomorphic encryption over the Integers with Shorter Public Key by Coron et.al
// http://publications.uni.lu/bitstream/10993/12396/1/441.pdf
// #define SEC_PARAM 42		// toy 
// #define SEC_PARAM 52		// small
// #define SEC_PARAM 62		// medium 
// #define SEC_PARAM 72		// large

#define SEC_PARAM 140
#define NR_TESTE 1000
// #define MULT_DEPTH 50

using namespace std;

vector<int> compute_random_circuit()
{
    vector<int> circuit(NR_TESTE);
    
    srand(time(NULL));
    for(int i=0; i<NR_TESTE; i++)
    {
        circuit[i] = rand() % 2 ;
    }
    
    return circuit;
}

void test_Flat_DGHV_simple()
{
	cout << "Testare Flat_DGHV_simple ...\n\n";
	int lambda = SEC_PARAM;

    // const char *source_file = "../Data/DGHVparams42.txt";
	// Flat_DGHV fdghv(source_file);
	Flat_DGHV fdghv(lambda);

	Mat_ZZ C1;
	Mat_ZZ C2;
	int m1;
	int m2;
    double wtime = 0.0;
    
    cout << "\nTestarea adancimii circuitului evaluat ...\n\n";
    
    ZZ c1, c2;  // ctxts of DGHV original scheme
    int max_depth = NR_TESTE;
    double mean_time = 0.0;
    vector<int> circuit = compute_random_circuit();

	for (int test = 0; test < 100; test++)
	{
        // vector<int> circuit = compute_random_circuit();
        max_depth = NR_TESTE;
             
        // simple DGHV;
        srand(time(NULL)); 
        m1 = rand() % 2;
        c1 = fdghv.encrypt_DGHV(m1);
        
        #if defined(_OPENMP)
                wtime = omp_get_wtime ( );       
        #endif

		for (int i = 0; i < NR_TESTE; i++)
		{
			m2 = rand() % 2;
			c2 = fdghv.encrypt_DGHV(m2);
            
            // if( circuit[i] == 1 )
            // {
                c1 = c1 * c2;
                m1 *= m2;
            /*}
            else
            {
                c1 = c1 + c2;
                m1 = ( m1 + m2 ) % 2;
            }*/
                 

			if ( fdghv.decrypt_DGHV(c1) != m1 )
			{
                max_depth = i;
				cout << "\tMax eval depth FDGHV i = " << i << endl;
				break;
			}

		}
        
         #if defined(_OPENMP)
		    wtime = omp_get_wtime ( ) - wtime;
        #endif
        mean_time += wtime;
        continue ;
        
        cout << "  TIME DGHV original = " << wtime << "\n";  
        
        
        m1 = 1;
		C1 = fdghv.encrypt(m1);

#if defined(_OPENMP)
        wtime = omp_get_wtime ( );
#endif
        int i = -1;
        bool passed = true;
		for (i = 0; i < NR_TESTE; i++)
		{
            // m2 = 1;
			m2 = rand() % 2;
		    C2 = fdghv.encrypt(m2);
            // C1 = fdghv.omp_hom_mult_opt(C1, C2);
            // C1 = fdghv.hom_mult_opt(C1, C2);
            m1 = m1 * m2;
            
            // C1 = fdghv.omp_hom_add(C1, C2);
            // m1 = ( m1 + m2 ) % 2;
             
            if( circuit[i] == 1 )
            {
                C1 = fdghv.omp_hom_mult(C1, C2);
                m1 *= m2;
                // cout <<  circuit[i] << " ";
            }
            else
            {
                C1 = fdghv.omp_hom_add(C1, C2);
                m1 = ( m1 + m2 ) % 2;
                // cout <<  circuit[i] << " ";
            }
            
			if ( fdghv.decrypt(C1) != m1 )
			{
				cout << " EVAL_FDGHV depth = ";
                cout << i << endl << endl;
                passed = false;
				break;
			}

		}
        
        #if defined(_OPENMP)
		    wtime = omp_get_wtime ( ) - wtime;
        #endif
        cout << "  TIME CU OMP = " << wtime << "\n\n";
        if( passed == true )
            cout << " test " << test << " - PASSED\n\n";
        mean_time += wtime;
        
        break;
  
		/*m1 = m2 = 1;
		C1 = C2 = fdghv.encrypt(m1);

#if defined(_OPENMP)
		wtime = omp_get_wtime ( );
#endif

		for (int i = 0; i < NR_TESTE; i++)
		{
			m2 = rand() % 2;
			C2 = fdghv.encrypt(m2);
            if( circuit[i] == 1 )
            {
                C1 = fdghv.hom_mult_opt(C1, C2);
                m1 *= m2;
            }
            else
            {
                C1 = fdghv.hom_add(C1, C2);
                m1 = ( m1 + m2 ) % 2;
            }
                 

			if ( fdghv.decrypt(C1) != m1 )
			{

				cout << "Max mult depth FDGHV i = " << i << endl;
				break;
			}

		}

#if defined(_OPENMP)
		wtime = omp_get_wtime ( ) - wtime;
#endif
        cout << "  TIME FARA OMP = " << wtime << "\n"; */

	}
    
    cout << "media = " << mean_time / 100 << endl;

	cout << "Final test Flat_DGHV_simple\n";

}

void test_FDGHV_max_depth()
{
    int lambda = SEC_PARAM;
    // const char *source_file = "../Data/DGHVparams42.txt";
	// Flat_DGHV fdghv(source_file);
	Flat_DGHV fdghv(lambda);
	Mat_ZZ C1;
	Mat_ZZ C2;
	int m1;
	int m2;
    double wtime = 0.0;
    vector<int> circuit = compute_random_circuit();

	m1 = 1;
	C1 = fdghv.encrypt(m1);

	cout << "Testing FDGHV scheme max depth ...\n";
	
    #if defined(_OPENMP)
        wtime = omp_get_wtime ( );

        LL i = 0;
        
        for (i = 0; i < NR_TESTE; i++)
        {
            m2 = rand() % 2;
            C2 = fdghv.encrypt(m2);
            
            if( circuit[i] == 1 )
            {
                C1 = fdghv.omp_hom_mult(C1, C2);
                m1 *= m2;
            }
            else
            {
                C1 = fdghv.omp_hom_add(C1, C2);
                m1 = ( m1 + m2 ) % 2;
            }
            
            if (fdghv.decrypt(C1) != m1 )
            {
                cout << "Max mult depth FDGHV = " << i << endl;
                break;
            }

        }
        
		wtime = omp_get_wtime ( ) - wtime;
    #endif
    cout << "  TIME CU OMP = " << wtime << "\n\n";
	cout << "Test pentru max_FDGHV_depth incheiat\n";
}

Flat_DGHV* HE_Context;
Mat_ZZ* ctxt_of_1;                  // encryption of 1 (constant)
int		t_bits;

void test_max_func()
{
    int lambda = SEC_PARAM;
    // const char *source_file = "../Data/DGHVparams42.txt";
	// Flat_DGHV fdghv(source_file);
	HE_Context = new Flat_DGHV(lambda);
    ctxt_of_1 = new Mat_ZZ( HE_Context->encrypt(1) );
    t_bits = 32;
    int max_val = (int)( pow(2,32) - 1 );
    
    vector<int> values;
    srand(time(NULL));
    for(int i=0; i<16; i++)
    {
        values.push_back( rand() % 1024);
    }
    
    vector<vector<Mat_ZZ> > vvct;
    
    HE_Integer hint((*HE_Context), t_bits);
    hint.encryptIntVector(vvct, values);
    
    clock_t tStart = clock();
    vector<Mat_ZZ> max = getmax(vvct, 0, vvct.size() );
    cout << "Timp : " << (double)(clock()-tStart)/CLOCKS_PER_SEC << endl;
    
    
    int plain_max = gmax(values, 0, 16);
    cout << "plain = " << plain_max << endl;
    cout << "hom = " << hint.decryptIntValue(max) << endl;
    
    // cleanup
    delete HE_Context;
    HE_Context = nullptr;
    delete ctxt_of_1;
    ctxt_of_1 = nullptr;
    t_bits = 0;
}

int main()
{
    test_max_func();
    
	// cout << " bitset pentru matricile GSW\n\n";
	
	// test_Flat_DGHV_simple();
    
    // test_FDGHV_max_depth();

	return 0;

}
