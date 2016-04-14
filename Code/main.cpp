# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <omp.h>

#include "utilities.h"
#include <time.h>
#include <assert.h>
#include "Flat_DGHV.h"
#include <fstream>
#include <NTL/ZZ.h>

// TODO : IMPLEMENTAREA SCHEMEI PE BAZA NTL
// DGHV securiy parameter
// dupa Fully homomorphic encryption over the Integers with Shorter Public Key by Coron et.al
// http://publications.uni.lu/bitstream/10993/12396/1/441.pdf
// valoare acceptata pentru x_0 si enc_0 reprezentati ca long
// #define SEC_PARAM 42		// toy 
// #define SEC_PARAM 52		// small
// #define SEC_PARAM 62		// medium 
// #define SEC_PARAM 72		// large

#define SEC_PARAM 100
#define NR_TESTE 400
#define MULT_DEPTH 50

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

	
	// m1 = rand() % 2;
	m1 = 1;
	C1 = fdghv.encrypt(m1);

	m2 = 1;
	C2 = fdghv.encrypt(m2);
    
    cout << "\nTestare multiplicari ...\n\n";
    
    vector<int> circuit = compute_random_circuit();

	for (int test = 0; test < 20; test++)
	{
        srand(time(NULL));
        
		m1 = m2 = 1;
		C1 = fdghv.encrypt(m1);
		C2 = fdghv.encrypt(m2);

#if defined(_OPENMP)
		wtime = omp_get_wtime ( );
#endif

		for (int i = 0; i < NR_TESTE; i++)
		{
			// m2 = rand() % 2;
			// C2 = fdghv.encrypt(m2);
            if( circuit[i] == 1 )
            {
                C1 = fdghv.hom_mult(C1, C2);
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
			else
			{
				// m1 = (m1 + m2) % 2;
				// m1 *= m2;
				// cout << "Succes\n";
			}

		}

#if defined(_OPENMP)
		wtime = omp_get_wtime ( ) - wtime;
#endif
        cout << "  TIME FARA OMP = " << wtime << "\n"; 
        
		m1 = 1;
		C1 = fdghv.encrypt(m1);
		m2 = 1;
		C2 = fdghv.encrypt(m2);

#if defined(_OPENMP)
        wtime = omp_get_wtime ( );
#endif

		for (int i = 0; i < NR_TESTE; i++)
		{
			// m2 = rand() % 2;
			// C2 = fdghv.encrypt(m2);

			// C1 = fdghv.omp_hom_mult(C1, C2);
            // C1 = fdghv.omp_hom_mult_opt(C1, C2);
            
            if( circuit[i] == 1 )
            {
                C1 = fdghv.omp_hom_mult_opt(C1, C2);
                m1 *= m2;
            }
            else
            {
                C1 = fdghv.omp_hom_add(C1, C2);
                m1 = ( m1 + m2 ) % 2;
            }
                 

			if (fdghv.decrypt(C1) != m1 )
			{
				cout << "Max mult depth FDGHV i = " << i << endl;
				break;
			}
			else
			{
				// m1 = (m1 + m2) % 2;
				// m1 *= m2;
				// cout << "Succes\n";
			}

		}
        
        #if defined(_OPENMP)
		    wtime = omp_get_wtime ( ) - wtime;
        #endif
        cout << "  TIME CU OMP = " << wtime << "\n\n";

		// continue;

	}

	cout << "Final test Flat_DGHV_simple\n";

}


int main()
{
	cout << " bitset pentru matricile GSW\n\n";
	
	test_Flat_DGHV_simple();

	return 0;

}
