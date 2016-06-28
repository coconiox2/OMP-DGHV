#pragma once

#include <iostream>
#include <vector>
#include <NTL/ZZ.h>

using namespace std;
using namespace NTL;

/*
@brief This class implements FHE over Integers
		using the public key compression technique from
		"Public Key Compression and Modulus Switching for FHE over Integers" by Coron et.al.
		http://www.iacr.org/archive/eurocrypt2012/72370441/72370441.pdf
*/
class CNT_FHE
{
	int lambda; // the security parameter
	int ro; // the no of bits of noise for the public key integers
	int ro_prim; // the no of bits of noise used at encryption
	int eta; // bit size of the secret key integer
	int gamma; // the no of bits of the chi_i's integers of the public key
	int tau; // the size of the public key vector
	int alpha; // size of b_i's used at encryption
    int beta;
    
    ZZ base;

public:
    CNT_FHE() {}

	CNT_FHE(int lambda, int ro, int ro_prim, int eta, int gamma, int tau, int alpha);

	~CNT_FHE() {}
    
    void set_base(ZZ b) { base = b; }
    
    void set_parameters(int lambda, int ro, int ro_prim, int eta, int gamma, int tau, int alpha);

	void generate_keys(ZZ &sk, vector<ZZ> &pk)const;
    
    void quadratic_key_generation(ZZ &sk, vector<ZZ> &pk)const;

	void encrypt(vector<ZZ> &pk, int message, ZZ &ctxt)const;
    
    void quadratic_encryption(vector<ZZ> &pk, int message, ZZ &ctxt)const;

	void decrypt(ZZ &sk, int &decrypted_message, ZZ &ctxt)const;
};

ZZ sample_r(int ro);


