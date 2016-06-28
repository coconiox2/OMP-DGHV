#pragma once
#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <vector>
#include "CNT_FHE.h"

using namespace std;
using namespace NTL;

typedef unsigned int uint32;
typedef vector<vector<uint32> > Mat_uint32;
typedef vector<uint32> Vec_uint32;

/*
@brief clasa care implementeaza FHE peste intregi
	folosind tehnicile de Flattening descrise pentru schema GSW
*/
class Flat_CNT
{
	ZZ w; 
	vector<ZZ> v;
	long l;

	ZZ sk_CNT;
	vector<ZZ> pk_CNT;
    
	vector<ZZ> C_prim;

    CNT_FHE cnt_fhe;
    
    void    compute_FDGHV_settings();

	Vec_uint32	bitdecomp(vector<ZZ> &C_i)const;		// BitDecomp(a) = {a_0, a_1, ..., a_n} , a = a_0 + 2*a_1+... +2^(n)*n

	vector<ZZ>	bitdecomp_1(Vec_uint32 &C_i)const;		// BitDecomp(a_0, a_1, ..., a_n) = a , cu a = a_0 + 2*a_1+... +2^(n)*n

	Mat_uint32	flatten(Mat_uint32 &C)const;			// Flatten(C)=BitDecomp( BitDecomp_1(C) % x_0 )

public:

	Flat_CNT(int lambda, ZZ baza = ZZ(2) );
    
    ~Flat_CNT() {}
    
	Mat_uint32	encrypt(int message)const;

	ZZ		decrypt(const Mat_uint32 &C)const;
    
    Mat_uint32	hom_mult_opt(Mat_uint32 &C1, Mat_uint32 &C2)const;
    
    vector<vector<ZZ> >    encrypt_ZZ(int m)const;
    int       decrypt_ZZ(vector<vector<ZZ> > C)const;

	/*Mat_uint32	hom_add(Mat_uint32 &C1, Mat_uint32 &C2)const;

	Mat_uint32	hom_mult(Mat_uint32 &C1, Mat_uint32 &C2)const;

	long get_l()const { return l; }

	Mat_uint32	hom_mult_opt(Mat_uint32 &C1, Mat_uint32 &C2)const;

    ZZ		encrypt_CNT(int message);

	ZZ		decrypt_CNT(ZZ &ctxt)const;
    
    ZZ      get_x_0()const { return pk_CNT[0]; }*/
};

