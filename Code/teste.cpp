#include "teste.h"
#include <fstream>
#include <ctime>

void test_max_mult_depth()
{
	ofstream out("max_depth.txt");

	double elapsed_time = 0.0;

	Mat_ZZ C1;
	Mat_ZZ C2;
	int m1;
	int m2;

	int decomp_base[] = { 32678, 128, 1024, 8192, 2};
	int security[] = { 112, 128, 140, 192 };

	int w_baza = -1;
	int lambda = -1;

	for (int i = 0; i < 5; i++)
	{
		w_baza = decomp_base[i];
		baza = new ZZ(w_baza);

		for (int j = 0; j < 3; j++)
		{
			lambda = security[j];
			
			Flat_DGHV he_engine(lambda, (*baza));

			cout << "SEC = " << lambda << " \t BAZA = " << w_baza << endl;
			out << "SEC = " << lambda << " \t BAZA = " << w_baza << endl;

			// test max mult depth
			m1 = w_baza - 1;
			C1 = he_engine.encrypt(m1);

			int mult_depth = 0;
		    clock_t begin2 = clock();
			do
			{
				C1 = he_engine.hom_mult_opt(C1, C1);
                m1 = ( m1 * m1 ) % w_baza;
			} while (he_engine.decrypt(C1) == m1);
            
            clock_t end2 = clock();
            elapsed_time = double(end2 - begin2) / CLOCKS_PER_SEC;
            
			cout << "Time Max Mult Depth = " << elapsed_time << endl;
			out << "Time Max Mult Depth = " << elapsed_time << endl;
			cout << "Max Mult Depth = " << mult_depth << endl << endl;
			out << "Max Mult Depth = " << mult_depth << endl << endl;
            
            cout << "***********************************************************" << endl << endl;
            out << "***********************************************************" << endl << endl;
		}
	}

	out.close();
}



void test_schema_FDGHV()
{
	ofstream out("masuratori.txt");

	double elapsed_time = 0.0;

	Mat_ZZ C1;
	Mat_ZZ C2;
	int m1;
	int m2;

	int decomp_base[] = { 2, 16, 128, 1024, 32678 };
	int security[] = { 192, 220, 140 };

	int w_baza = -1;
	int lambda = -1;

	for (int i = 0; i < 5; i++)
	{
		w_baza = decomp_base[i];
		baza = new ZZ(w_baza);

		for (int j = 0; j < 3; j++)
		{
			lambda = security[j];
			
			Flat_DGHV he_engine(lambda, (*baza));

			cout << "SEC = " << lambda << " \t BAZA = " << w_baza << endl;
			out << "SEC = " << lambda << " \t BAZA = " << w_baza << endl;

			srand(time(NULL));

			m1 = rand() % w_baza;
			m2 = rand() % w_baza;

			C1 = he_engine.encrypt(m1);
			C2 = he_engine.encrypt(m2);

			bool circuited_completed = true;
			int circuit_size = 100;

			for (int l = 0; l < 3; l++)
			{
				int mults = 0;
				int adds = 0;
				vector<int> circuit = compute_random_circuit(circuit_size);

				clock_t begin = clock();
				for (int k = 0; k < circuit_size; k++)
				{
					if (circuit[k] == 0)
					{
						C1 = he_engine.hom_add(C1, C2);
						m1 = (m1 + m2) % w_baza;
						adds++;
					}
					else
					{
						C1 = he_engine.hom_mult_opt(C1, C2);
						m1 = (m1*m2) % w_baza;
						mults++;
					}

					if (he_engine.decrypt(C1) != m1)
					{
						cout << "Eval depth = " << k << endl;
						out << "Eval depth = " << k << endl;
						circuited_completed = false;
						break;
					}
					m2 = rand() % w_baza;
					C2 = he_engine.encrypt(m2);
				}
                clock_t end = clock();
                elapsed_time = double(end - begin) / CLOCKS_PER_SEC;
				cout << "Time = " << elapsed_time << " s" << endl;
				out << "Time = " << elapsed_time << " s" << endl;

				cout << "Mults = " << mults << " Adds = " << adds << endl;
				out << "Mults = " << mults << " Adds = " << adds << endl;

				if (circuited_completed == true)
				{
					cout << "Circuit completed depth = " << circuit_size << endl << endl;
					out << "Circuit completed depth = " << circuit_size << endl << endl;
				}
				else
				{
					break;
				}
				circuit_size *= 10;
			}

			// test max mult depth
            /* cout<< "Se testeaza max mult depth ..." << endl;
			srand(time(NULL));
			m1 = 1;
			C1 = he_engine.encrypt(m1);

			int mult_depth = 0;
		    clock_t begin2 = clock();
			do
			{
				C1 = he_engine.hom_mult_opt(C1, C1);
			} while (he_engine.decrypt(C1) == 1);
            clock_t end2 = clock();
            elapsed_time = double(end2 - begin2) / CLOCKS_PER_SEC;
			cout << "Time Max Mult Depth = " << elapsed_time << endl;
			out << "Time Max Mult Depth = " << elapsed_time << endl;
			cout << "Max Mult Depth = " << mult_depth << endl << endl;
			out << "Max Mult Depth = " << mult_depth << endl << endl;*/
            
            cout << "***********************************************************" << endl << endl;
            out << "***********************************************************" << endl << endl;
		}
	}

	out.close();
}


vector<int> compute_random_circuit(int size)
{
    vector<int> circuit(size);
    
    srand(time(NULL));
    for(int i=0; i<size; i++)
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

void test_max_no_encryptions()
{
    int w_baza = 2;
	int lambda = 112;
	baza = new ZZ(w_baza);
    
    Flat_DGHV he_engine(lambda, (*baza));
    
    vector<Mat_ZZ> encs;
    int i = 0;
    
    try
    {
        while(true)
        {
            Mat_ZZ C = he_engine.encrypt(i%2);
            encs.push_back(C);
            i++;
            cout<< i << " ";
        }
    }
    catch (std::bad_alloc& ba)
    {
        cerr << "bad_alloc caught: " << ba.what() << endl;
        delete baza;
        baza = nullptr;
    }
    catch(...)
    {
        cout << "Exceptie ^_^ ..." << endl;
        delete baza;
        baza = nullptr;
    }
    
    
}

// int rep_baza[] = { 1, 2, 4, 5, 7, 9, 10, 15};    
// Nu putem beneficia de avantajul bazelor de descompunere mai mari 
// intrucat avem nevoie de operatiile binare AND si XOR pentru 
// a construi circuitul de evaluare al functiei getmax

// valorile implicite : lambda = 192, no_bits = 8, N = 16
void test_max_func(int lambda, int no_bits, int N)
{
    ofstream out;
    out.open ("getmax.txt", std::ofstream::out | std::ofstream::app);
    
    /**** parametrii schemei pentru testare ****/
    t_bits = no_bits;
    int w_baza = 2;
                    
    cout<< endl;
    cout<< "security_bits = " << lambda << endl;
    cout<< "int_size = " << t_bits << endl;
    cout<< "N = " << N << endl; 
    cout<< "baza = " << w_baza << endl;
    
    double average_time = 0.0;
                    
    baza = new ZZ(w_baza);
    cout << "Se creeaza contextul homomorfic ... " << endl;
    HE_Context = new Flat_DGHV(lambda, (*baza));
    cout << "HE_Context creat." << endl;
    ctxt_of_1 = new Mat_ZZ( HE_Context->encrypt(1) );
    cout << "Se creaza HE_Integer ... " << endl;
    HE_Integer hint((*HE_Context), t_bits);
    cout << "HE_Integer creat." << endl;
    vector<vector<Mat_ZZ> > vvct;

    int max_val = (int)( pow(2, t_bits) - 1 );
    vector<int> values(N);
                
    srand(time(NULL));
    for(int i=0; i<N; i++)
    {
        values[i] = rand() % max_val;
    }

    try
    {
        cout << "Encrypting int vector ..." << endl;
        hint.encryptIntVector(vvct, values);
        cout << "vvct.size = " << vvct.size() << endl;
    }    
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc  encryptIntVector caught: " << ba.what() << '\n';
    }
    catch(...)
    {
        cout << "Exceptie straina ..." << endl;
    }
                    
    vector<Mat_ZZ> max; // = getmax(vvct, 0, vvct.size() );
    clock_t tStart = clock();
    try
    {
        cout << "A intrat in getmax." << endl;
        max = getmax(vvct, 0, vvct.size() );
        cout << "A iesit din getmax." << endl;
    }    
    catch (std::bad_alloc& ba)
    {
        std::cerr << "bad_alloc getmax caught: " << ba.what() << '\n';
    }
    cout << "Timp homomorphic : " << (double)(clock()-tStart)/CLOCKS_PER_SEC << endl;

    out<< "security_bits = " << lambda << endl;
    out<< "int_size = " << t_bits << endl;
    out<< "N = " << N << endl; 
    out<< "baza = " << w_baza << endl;
    out << "Timp homomorphic : " << (double)(clock()-tStart)/CLOCKS_PER_SEC << endl;
    average_time += (double)(clock()-tStart)/CLOCKS_PER_SEC;
                    
    int plain_max = gmax(values, 0, values.size());
    cout << "plain_max = " << plain_max << endl;
    cout << "hom_max = " << hint.decryptIntValue(max) << endl << endl;
    cout << "**********************************************************************" << endl << endl;
                    
    // writing to file
    out << "plain_max = " << plain_max << endl;
    out << "hom_max = " << hint.decryptIntValue(max) << endl << endl;
    out << "**********************************************************************" << endl << endl;
    
    // cleanup
    delete HE_Context;
    HE_Context = nullptr;
    delete ctxt_of_1;
    ctxt_of_1 = nullptr;
    delete baza;
    baza = nullptr;
    
    t_bits = 0;
    out.close();
}


void cmp_HElib()
{
    ofstream out;
    out.open ("cmp_HElib.txt", std::ofstream::out | std::ofstream::app);
    
    /**** parametrii testare ****/
    int w_baza = 2;
    baza = new ZZ(w_baza);
    
    int security_bits[] = { 48, 106, 140};
    t_bits = 8;
    int N = 16;
    
    // int rep_baza[] = { 1, 2, 4, 5, 7, 9, 10, 15};    
    // Nu putem beneficia de avantajul bazelor de descompunere mai mari 
    // intrucat avem nevoie de operatiile binare AND si XOR pentru 
    // a construi circuitul de evaluare al functiei getmax
    /**** parametri testare ****/
    cout<< "int_size = " << t_bits << endl;
    cout<< "N = " << N << endl; 
    cout<< "baza = " << w_baza << endl;
    
    out<< "int_size = " << t_bits << endl;
    out<< "N = " << N << endl; 
    out<< "baza = " << w_baza << endl;
    
    
    for(int sb=0; sb<3; sb++)
    {
        int lambda = security_bits[sb];

        double average_time = 0.0;       
                    
        HE_Context = new Flat_DGHV(lambda, (*baza));
        ctxt_of_1 = new Mat_ZZ( HE_Context->encrypt(1) );
        HE_Integer hint((*HE_Context), t_bits);
        vector<vector<Mat_ZZ> > vvct;

        int max_val = (int)( pow(2, t_bits) - 1 );
        vector<int> values(N);
                
        srand(time(NULL));
        for(int i=0; i<N; i++)
        {
             values[i] = rand() % max_val;
        }

        try
        {
            hint.encryptIntVector(vvct, values);
        }    
        catch (std::bad_alloc& ba)
        {
             std::cerr << "bad_alloc  encryptIntVector caught: " << ba.what() << '\n';
        }
                    
        clock_t tStart = clock();
        vector<Mat_ZZ> max; // = getmax(vvct, 0, vvct.size() );
        try
        {
            max = getmax(vvct, 0, vvct.size() );
        }    
        catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc getmax caught: " << ba.what() << '\n';
        }
        cout << "Timp homomorphic : " << (double)(clock()-tStart)/CLOCKS_PER_SEC << endl;
        out << "Timp homomorphic : " << (double)(clock()-tStart)/CLOCKS_PER_SEC << endl;
        average_time += (double)(clock()-tStart)/CLOCKS_PER_SEC;
                    
        int plain_max;
        tStart = clock();
        plain_max = gmax(values, 0, values.size());
        cout << "Timp plain : " << (double)(clock()-tStart) << endl;
        cout << "plain_max = " << plain_max << endl;
        cout << "hom_max = " << hint.decryptIntValue(max) << endl << endl;
        cout << "**********************************************************************" << endl << endl;
                    
        // writing to file
        out << "Timp plain : " << (double)(clock()-tStart)/1000000.0F << endl;
        out << "plain_max = " << plain_max << endl;
        out << "hom_max = " << hint.decryptIntValue(max) << endl << endl;
        out << "**********************************************************************" << endl << endl;
                    
        // cleanup
        delete HE_Context;
        HE_Context = nullptr;
        delete ctxt_of_1;
        ctxt_of_1 = nullptr;
    }
    
    delete baza;
    baza = nullptr;
    t_bits = 0;
    out.close();
  
}