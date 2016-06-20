#include "HE_Utils.h"
#include <assert.h>
#include <fstream>

extern Flat_DGHV* HE_Context;
extern Mat_ZZ* ctxt_of_1;                  // encryption of 1 (constant)
extern int		t_bits;

// #define _OMP_USAGE_


/*************************************************************************************/
Mat_ZZ compute_z(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
        Mat_ZZ ret;
        
        #ifdef _OMP_USAGE_
            ret = HE_Context->omp_hom_add( (*ctxt_of_1), ct_x[i]);
            ret = HE_Context->omp_hom_add(ret, ct_y[i]);
        #else
            ret = HE_Context->hom_add( (*ctxt_of_1), ct_x[i]);
            ret = HE_Context->hom_add(ret, ct_y[i]);
        #endif
        
		return ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

    Mat_ZZ ret = compute_z(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct = compute_z(i, l, ct_x, ct_y);
    
	ret = HE_Context->hom_mult_opt(ret, ct);                // ret *= ct;	

	return ret;
}

/*************************************************************************************/
Mat_ZZ compute_t(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
        Mat_ZZ ret = ct_x[i];
		ret = HE_Context->hom_mult_opt(ret, ct_y[i]);  // ret *= ct_y[i];
        #ifdef _OMP_USAGE_
		    ret = HE_Context->omp_hom_add(ret, ct_x[i]);       // ret += ct_x[i];
        #else
            ret = HE_Context->hom_add(ret, ct_x[i]);
        #endif
        
		return ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

    Mat_ZZ ret = compute_t(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct_z = compute_z(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct_t = compute_t(i, l, ct_x, ct_y);

	ct_z = HE_Context->hom_mult_opt(ct_z, ct_t);  // ct_z *= ct_t;
    #ifdef _OMP_USAGE_
        ret = HE_Context->omp_hom_add(ret, ct_z);
    #else
	    ret = HE_Context->hom_add(ret, ct_z);		// ret += ct_z;
    #endif

	return ret;
}

/*************************************************************************************/
Mat_ZZ compute_s(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		Mat_ZZ ret = ct_x[i];
		ret = HE_Context->hom_mult_opt(ret, ct_y[i]);		// ret *= ct_y[i];
        #ifdef _OMP_USAGE_
            ret = HE_Context->omp_hom_add(ret, ct_y[i]);			// ret += ct_y[i];
		    ret = HE_Context->omp_hom_add(ret, (*ctxt_of_1));
        #else
		    ret = HE_Context->hom_add(ret, ct_y[i]);			// ret += ct_y[i];
		    ret = HE_Context->hom_add(ret, (*ctxt_of_1));		// ret += *Mat_ZZ_of_1;
        #endif
		return ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	Mat_ZZ ret = compute_t(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct_z = compute_z(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct_s = compute_s(i, l, ct_x, ct_y);

	ct_z = HE_Context->hom_mult_opt(ct_z, ct_s);	// ct_z *= ct_s;
    #ifdef _OMP_USAGE_
        ret = HE_Context->omp_hom_add( ret, ct_z);
    #else
	    ret = HE_Context->hom_add(ret, ct_z);		    // ret += ct_z;
    #endif

	return ret;
}

/*************************************************************************************/
vector<Mat_ZZ> select(Mat_ZZ& c, vector<Mat_ZZ>& a, vector<Mat_ZZ>& b)
{
	vector<Mat_ZZ> ret;

	// HE_Integer heInt(*parms, t_bits);
	// cout << "a = " << heInt.decryptIntValue(a, *secret_key) << endl;
	// cout << "b = " << heInt.decryptIntValue(b, *secret_key) << endl;

	vector<Mat_ZZ> vt1, vt2;
	for (int i = 0; i < a.size(); i++)
	{
		vt1.push_back(Mat_ZZ(c));
		vt2.push_back(Mat_ZZ(c));
	}

	for (int i = 0; i < a.size(); i++)
	{
		vt1[i] = HE_Context->hom_mult_opt(vt1[i], a[i]);     // vt1[i] *= a[i];
        
        #ifdef _OMP_USAGE_
            vt2[i] = HE_Context->omp_hom_add(vt2[i], *ctxt_of_1);
        #else
		    vt2[i] = HE_Context->hom_add(vt2[i], *ctxt_of_1);    // vt2[i] += ENC(1);
        #endif
        
		vt2[i] = HE_Context->hom_mult_opt(vt2[i], b[i]);     // vt2[i] *= b[i];
            
        #ifdef _OMP_USAGE_
            vt1[i] = HE_Context->omp_hom_add(vt1[i], vt2[i]);
        #else
		    vt1[i] = HE_Context->hom_add(vt1[i], vt2[i]);        // vt1[i] += vt2[i];
        #endif
		ret.push_back(vt1[i]);
	}

	return ret;
}

/*************************************************************************************/

vector<Mat_ZZ> getmax(vector<vector<Mat_ZZ> >& vvct, int start, int n)
{
	assert(n >= 1);

	if (n == 1) return vvct[start];

	vector<Mat_ZZ> ct_max_1 = getmax(vvct, start, n / 2);
	vector<Mat_ZZ> ct_max_2 = getmax(vvct, start + n / 2, n % 2 == 0 ? n / 2 : n / 2 + 1);

	Mat_ZZ ct_t = compute_t(0, t_bits, ct_max_1, ct_max_2);
	
	vector<Mat_ZZ> ct_max = select(ct_t, ct_max_1, ct_max_2);

	return ct_max;

}

int gmax(vector<int>& v, int start, int n)
{
	assert(n >= 1);
	if (n == 1) return v[start];

	int max_1 = gmax(v, start, n / 2);
	int max_2 = gmax(v, start + n / 2, n % 2 == 0 ? n / 2 : n / 2 + 1);

	return max_1 > max_2 ? max_1 : max_2;
}
