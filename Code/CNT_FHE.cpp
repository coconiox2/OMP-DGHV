#include "CNT_FHE.h"

CNT_FHE::CNT_FHE(int l, int r, int rp, int e, int g, int t, int a)
{
	set_parameters(l, r, rp, e, g, t, a);
}

void CNT_FHE::set_parameters(int l, int r, int rp, int e, int g, int t, int a)
{
    lambda = l;
    ro = r;
    ro_prim = rp;
    eta = e;
    gamma = g;
    tau = t;
    alpha = a;
    
    base = ZZ(2);
    
    cout << "Parametri Schema CNT:" << endl;
	cout << "lambda = " << lambda << endl;
	cout << "ro = " << ro << endl;
	cout << "ro_prim = " << ro_prim << endl;
	cout << "eta = " << eta << endl;
	cout << "gamma = " << gamma << endl;
	cout << "tau = " << tau << endl;
	cout << "alpha = " << alpha << endl << endl;
}

void CNT_FHE::quadratic_key_generation(ZZ &p, vector<ZZ> &pk)const
{
    int beta = tau;
    
    RandomBits(p, eta);
    
    pk = vector<ZZ>(1+2*beta);
    
    ZZ q_0; // q_0 from [0, 2^gamma/p]
	ZZ max = power(ZZ(2), gamma) / p;
    do
    {
	    RandomBnd(q_0, max);
    }while( ( q_0 % ZZ(2) == 1 ) && (q_0 % base == 1) ) ;
    
	pk[0] = p*q_0; // x_0 = p * q_0
    
    ZZ q_i;
    ZZ r_i;
    // generate x_i,b for 1<=i<=beta and b from {0, 1}
    for(int i=0; i<2*beta; i++)
    {
        RandomBnd(q_i, q_0);
        r_i = sample_r(ro);
        pk[1+i] = p*q_i + r_i;
    }
}

void CNT_FHE::quadratic_encryption(vector<ZZ> &pk, int message, ZZ &ctxt)const
{
    cout << "Criptarea nu este implementate corect." << endl;
    int beta = tau;
    cout << "beta = " << beta << endl;
    
    vector<ZZ> b(beta*beta);
	ZZ two_alpha = power(ZZ(2), alpha);
    
	for (int i = 0; i < beta*beta; i++)
	{
		RandomBnd(b[i], two_alpha);
	}
    
    ZZ sum(0);
    
    for(int i=0; i<beta; i++)
    {
        for(int j=0; j<beta; j++)
        {
            // sum += b_i,j * x_i,0 * x_j,1
            sum += b[i*beta+j] * pk[1+i] * pk[1+beta+j];
            break;
        }
        sum = sum * 44;
        break;
        // cout << "checkpoint i = " << i << endl;
    }
    
    // cout<<"for2\n";
    
    ZZ r = sample_r(ro_prim);
    
    ctxt = ( ZZ(message) + base*r + base*sum ) % pk[0];
}

void CNT_FHE::generate_keys(ZZ &p, vector<ZZ> &pk)const
{
	RandomBits(p, eta);

	pk = vector<ZZ>(tau+2);

	ZZ seed;
	RandomBits(seed, 32);
	pk[0] = seed;

	ZZ q_0; // q_0 from [0, 2^gamma/p]
	ZZ max = power(ZZ(2), gamma) / p;
	RandomBnd(q_0, max);
	pk[1] = p*q_0; // x_0 = p * q_0


	vector<ZZ> zeta_i(tau);
	ZZ max_zeta = power(ZZ(2), lambda + eta) / p;
	vector<ZZ> r_i(tau);
	for (int i = 0; i < tau; i++)
	{
		RandomBnd(zeta_i[i], max_zeta);
		r_i[i] = sample_r(ro);
	}

	ZZ chi_i;
	ZZ delta_i;
	ZZ x_i;
    ZZ max_chi = power(ZZ(2), gamma);

	SetSeed(seed);
	for (int i = 0; i < tau; i++)
	{
		// RandomBits(chi_i, gamma);
        RandomBnd(chi_i, max_chi);

		delta_i = chi_i % p + zeta_i[i]*p - r_i[i];

		pk[i + 2] = delta_i;
	}
}

void CNT_FHE::encrypt(vector<ZZ> &pk, int message, ZZ &ctxt)const
{
	vector<ZZ> b(tau);
	ZZ two_alpha = power(ZZ(2), alpha);
    
	for (int i = 0; i < tau; i++)
	{
		RandomBnd(b[i], two_alpha);
	}

	ZZ r = sample_r(ro_prim);
    
	ZZ sum(0);
	ZZ chi_i;
    ZZ max_chi = power(ZZ(2), gamma);

    SetSeed(pk[0]);
    for (int i = 0; i < tau; i++)
    {
        // RandomBits(chi_i, gamma);
        RandomBnd(chi_i, max_chi);
        sum += b[i]*(chi_i - pk[2 + i]);
    }

	ctxt = (ZZ(message) + 2 * r + 2 * sum) % pk[1];
}

void CNT_FHE::decrypt(ZZ &sk, int &decrypted_message, ZZ &ctxt)const
{
	conv( decrypted_message, ctxt % sk % base);
}

ZZ sample_r(int ro)
{
	ZZ r(0);

	// RandomBits(r, 2*ro);
    ZZ max = power(ZZ(2), 2*ro);
    RandomBnd(r, max);

	ZZ two_ro = power(ZZ(2), ro);

	r = r - two_ro;

	return r;
}

