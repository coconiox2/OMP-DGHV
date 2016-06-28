#include "Flat_CNT.h"
#include "Params.h"
#include "utilities.h"
#include <assert.h>

Flat_CNT::Flat_CNT(int lambda, ZZ baza)
{
	// w = ZZ(32 * 1024);
    w = ZZ(2);
    
    int l = 16;
	int r = l;
	int e = (int)pow(l, 2);
	int g = (int)pow(l, 2);
	int b = l*l;
	int a = l; // alpha * beta^2 >= gama + w(log lambda)
    int rp = 4*l;
    
    
    /*int l = 42;
	int r = 16;
	int e = 336;
	int g = 45000;
	int b = 12;
	int a = 1000; // alpha * beta^2 >= gama + w(log lambda)
    int rp = 4*l;*/
    
    cnt_fhe.set_parameters(l, r, rp, e, g, b, a);
    // cnt_fhe.set_base(w);
    
    cout << "Flat_CNT::Flat_CNT se genereaza cheile CNT_FHE ..." << endl;
    cnt_fhe.quadratic_key_generation(sk_CNT, pk_CNT);
    cout << "Chei generate." << endl;
    
    compute_FDGHV_settings();
}

void Flat_CNT::compute_FDGHV_settings()
{
	l = 0;					// l = log x_0 + 1
	ZZ x_0 = pk_CNT[0];
    cout << "Se calculeaza l ..." << endl;
	while (x_0 != 0)
	{
		l++;
		x_0 = x_0 / w;
	}
    
    cout << " l = " << l << endl;
    
	assert(l != 0);
	l += 1;
    cout << "l calculat." << endl;

	v = vector<ZZ>(l);		// v = Powersof2(1);
	ZZ two = w;
	ZZ pow_of_w(1);
    cout << "se creeaza vectorul v ..." << endl;
	for (int i = 0; i < l; i++)
	{
        v[i] = pow_of_w;
		pow_of_w *= w;
	}
    cout << "v creat." << endl;

	C_prim = vector<ZZ>(l);
    cout << "Generare C_0 incorecta." << endl;
    ZZ enc_of_0;
    cout << "Se cripteaza 0 ..." << endl;
    cnt_fhe.quadratic_encryption(pk_CNT, 0, enc_of_0);
    cout << "Zero criptat." << endl;
	for (int i = 0; i < l; i++)
	{
        // cnt_fhe.quadratic_encryption(pk_CNT, 0, C_prim[i]);
		C_prim[i] = enc_of_0; 
	}
}

/*ZZ	Flat_CNT::encrypt_CNT(int message)
{
	ZZ ctxt;
    cnt_fhe.quadratic_encryption(pk_CNT, message, ctxt);
    return ctxt;
}

ZZ	Flat_CNT::decrypt_CNT(ZZ &ctxt)const
{
	return ctxt % sk_CNT % w ;
}*/

Vec_uint32 Flat_CNT::bitdecomp(vector<ZZ> &C_i)const
{
	long length = C_i.size() * l;
    Vec_uint32 C_decomp(length, 0);

	ZZ elem;
	for (int i = 0, j=-1; i < length; i++)
	{
		if (i % l == 0)
		{
			j++;
			elem = C_i[j];
		}

        conv(C_decomp[i], elem % w);
		elem = elem / w;
	}
    
	return C_decomp;
}

vector<ZZ> Flat_CNT::bitdecomp_1(Vec_uint32 &C_i)const
{
	long length = C_i.size() / l;
    vector<ZZ> C_decomp_1(length);

	ZZ pow_of_w(1);
	for (LL i = 0, j = -1; i < C_i.size(); i++)
	{
		if (i % l == 0)
		{
			pow_of_w = 1;
			j++;
            C_decomp_1[i] = ZZ(0);		//C_decomp_1[j] = 0;
		}

		C_decomp_1[j] +=  ZZ(C_i[i]) * pow_of_w;
		pow_of_w *= w;
	}

	return C_decomp_1;
}

Mat_uint32 Flat_CNT::flatten(Mat_uint32 &C)const
{
	Mat_uint32 Flat_C(l);
	vector<ZZ> bd_1(1);
    
    LL i;
    ZZ x_0 = pk_CNT[0];
    int C_size = C.size();

    for (i = 0; i < C_size; i++)
    {
        bd_1 = bitdecomp_1(C[i]);
        bd_1[0] = bd_1[0] % x_0;
        Flat_C[i] = bitdecomp(bd_1);
    }

	return Flat_C;
}

Mat_uint32 Flat_CNT::encrypt(int message)const
{
    
    cout << "Criptare implementata incorect." << endl;
    
	Mat_uint32 C(l);
    
    vector<ZZ> linie(1, C_prim[0]);
	C[0] = bitdecomp(linie);
    int i;
    int dim = l;
#pragma omp parallel for \
    default(none) shared(C, message, dim) private(i)
	for (i = 1; i < dim; i++)
	{
        C[i] = C[0];
        C[i][i] += message;
	}
    C[0][0] += message;
    
    cout << "Intra in flatten." << endl;
    
	flatten(C);

	return C;
}

ZZ Flat_CNT::decrypt(const Mat_uint32 &C)const
{
	ZZ message(0);
	for (int i = 0; i < l; i++)
	{
		message += ZZ(C[0][i]) * v[i];
	}
	return message % sk_CNT % w ;
}

Mat_uint32 Flat_CNT::hom_mult_opt(Mat_uint32 &C1, Mat_uint32 &C2)const
{
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_uint32 C_mult(l);
	vector<vector<ZZ> > C1_prim(l);
    
	for (int i = 0; i < l; i++)
	{
        C1_prim[i] = bitdecomp_1(C1[i]);
    }

	for (int i = 0; i < l; i++)
	{
		ZZ elem(0);
		
		for (int j = 0; j < l; j++)
		{
			elem += ZZ(C2[i][j]) * C1_prim[j][0] ;
		}

		vector<ZZ> linie(1, elem);
        C_mult[i] = bitdecomp(linie);
	}

	return flatten(C_mult);
}

/*Mat_uint32 Flat_CNT::hom_add(Mat_uint32 &C1, Mat_uint32 &C2)const
{
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_uint32 C_add(l);
	
	for (int i = 0; i < l; i++)
	{
		for (int j = 0; j < l; j++)
		{
            C_add[i][j] = C1[i][j] + C2[i][j];
		}
	}

	return flatten(C_add);
}

Mat_uint32 Flat_CNT::hom_mult(Mat_uint32 &C1, Mat_uint32 &C2)const
{
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_uint32 C_mult(l);

	ZZ z;
	for (int i = 0; i < l; i++)
	{
        C_mult[i] = vector<ZZ>(l);
		
		for (int j = 0; j < l; j++)
		{
            C_mult[i][j] = ZZ(0);
			for (int k = 0; k < l; k++)
			{
				C_mult[i][j] += C1[i][k] * C2[k][j];
			}
		}
	}

	return flatten(C_mult);
}

*/


/*****************				OMP						*******
Mat_uint32 Flat_CNT::omp_hom_add(Mat_uint32 &C1, Mat_uint32 &C2)const
{
#if defined(_OPENMP)

	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_uint32 C_add(l);
	LL i, j;
    
#pragma omp parallel for collapse(2) \
    default(none) shared(C_add, C1, C2) private(i)
	for (i = 0; i < l; i++)
	{
		for (j = 0; j < l; j++)
		{
			C_add[i].push_back(C1[i][j] + C2[i][j]);
		}
    }
    
    return flatten(C_add);
    
#else
	return hom_add(C1, C2);
#endif

}

Mat_uint32 Flat_CNT::omp_hom_mult(Mat_uint32 &C1, Mat_uint32 &C2)const
{
#if defined(_OPENMP)
	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_uint32 C_mult(l);

	LL i, j, k;
	ZZ elem;

#pragma omp parallel for collapse(2) \
	default(none) shared(C_mult, C1, C2) private(i, j, k, elem) 
		for (i = 0; i<l; i++)
		{
            //#pragma omp parallel for collapse(2) \
	            default(none) shared(C_mult, C1, C2) \
                private(i, j, k, elem) 
			for (j = 0; j<l; j++)
			{
				elem = 0;
				for (k = 0; k<l; k++)
				{
					elem += C1[i][k] * C2[k][j];
                    //C_mult[i][j] += C1[i][k] * C2[k][j];
                }	
                C_mult[i].push_back(elem);			
			}
		}

	return flatten(C_mult);
#else
	return hom_mult(C1, C2);
#endif

	
}

/*
metoda abandonata deoarece nu este mai buna decat omp_hom_mult
pentru o paralelizare mai buna trebuie modificat algoritmul
de multiplicare combinat cu Flattening-ul
*/
/*Mat_uint32 Flat_CNT::omp_hom_mult_opt(Mat_uint32 &C1, Mat_uint32 &C2)const
{
#if defined(_OPENMP)

	assert(C1.size() == C2.size());
	assert(C1.size() == l);
	assert(C1[0].size() == l);

	Mat_uint32 C_mult(l);
    Mat_uint32 C1_prim(l);
	LL i, j, k;
    
	for (i = 0; i < l; i++)
	{
		C1_prim[i] = bitdecomp_1(C1[i]);
        C_mult[i].push_back( ZZ(0) );
	}
    
    //#pragma omp parallel for  \
        default(none) shared(C_mult, C1_prim, C2) private(i,j)
    for (i = 0; i < l; i++)
	{
		for (j = 0; j < l; j++)
		{
            C_mult[i][0] += C2[i][j] * C1_prim[j][0];
		}
	//}
    
    //for (i = 0; i < l; i++)
	//{
        C_mult[i] = bitdecomp(C_mult[i]);
    }

	return C_mult;

#else
	return hom_mult_opt(C1, C2);
#endif
}

Mat_uint32 Flat_CNT::omp_encrypt(int message)const
{
#if defined(_OPENMP)
	Mat_uint32 C(l);
	LL i;
	Vec_uint32 linie;

#pragma omp parallel for shared(C) private(i, linie)
	for (i = 0; i < l; i++)
	{
        linie = C_prim[i];
		C[i] = bitdecomp(linie);
		C[i][i] += message;
	}

	return flatten(C);

#else
	return encrypt(message);
#endif
}

int Flat_CNT::omp_decrypt(Mat_uint32 &C)const
{
#if defined(_OPENMP)

	ZZ message(0);
	LL i;
#pragma omp parallel for shared(message) private(i)
	for (i = 0; i < l; i++)
	{
		message += C[0][i] * v[i];
	}

	// return decrypt_DGHV(message);
    
    const ZZ miu = decrypt_DGHV(message);
    int val;
    conv(val, miu);
    return val;

#else
	return decrypt(C);
#endif
	
}


/*Vec_uint32 Flat_CNT::omp_bitdecomp(Vec_uint32 &C_i)const
{
    
#if defined(_OPENMP)
	long length = C_i.size() * l;
	Vec_uint32 C_decomp(length);
	ZZ elem;
	LL i;
    LL j = -1;
        
#pragma omp parallel for shared(C_decomp, C_i, j) private(elem, i)
	for (i = 0; i < length; i++)
	{
		
        #pragma omp critical
        {
            if (i % l == 0)
            {
                j++;
                elem = C_i[j];
            }
            C_decomp[i] = ZZ(elem % 2);
        }
		
		elem = elem / 2;
	}
    
	return C_decomp;
#else
	return bitdecomp(C_i);
#endif

}

Vec_uint32 Flat_CNT::omp_bitdecomp_1(Vec_uint32 &C_i)const
{
#if defined(_OPENMP)
	long length = C_i.size() / l;
	Vec_uint32 C_decomp_1(length);
	ZZ pow_of_two(1);
	ZZ two(2);
	UL i, j;
#pragma omp parallel for shared(C_decomp_1) private(i,j)
	for (i = 0, j = -1; i < C_i.size(); i++)
	{
		if (i % l == 0)
		{
			pow_of_two = 1;
			j++;
			C_decomp_1[i] = ZZ(0);		//C_decomp_1[j] = 0;
		}

		C_decomp_1[j] += C_i[i] * pow_of_two;
		pow_of_two *= two;
	}

	return C_decomp_1;
#else
	return omp_bitdecomp_1(C_i);
#endif

}*/