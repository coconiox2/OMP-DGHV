#include "HE_Utils.h"
// #include "HE_Integer.h"
#include <assert.h>
#include <fstream>

extern Flat_DGHV* HE_Context;
extern Mat_ZZ* ctxt_of_1;                  // encryption of 1 (constant)
extern int		t_bits;


/*************************************************************************************/
Mat_ZZ& compute_z(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		// ret = ctxt_of_1 + ct_x[i] + ct_y[i];
		Mat_ZZ* ret = new Mat_ZZ();
        *ret = HE_Context->hom_add( (*ctxt_of_1), ct_x[i]);
		*ret = HE_Context->hom_add(*ret, ct_y[i]);
		return *ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	Mat_ZZ* ret = new Mat_ZZ( compute_z(i + l, j - l, ct_x, ct_y) );
	Mat_ZZ ct = compute_z(i, l, ct_x, ct_y);

	*ret = HE_Context->hom_mult_opt(*ret, ct);                // ret *= ct;	
	return *ret;
}

/*************************************************************************************/
Mat_ZZ& compute_t(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		Mat_ZZ *ret = new Mat_ZZ(ct_x[i]);
		*ret = HE_Context->hom_mult_opt(*ret, ct_y[i]);  // ret *= ct_y[i];
		*ret = HE_Context->hom_add(*ret, ct_x[i]);       // ret += ct_x[i];
		return *ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	Mat_ZZ *ret = new Mat_ZZ( compute_t(i + l, j - l, ct_x, ct_y) );
	Mat_ZZ ct_z = compute_z(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct_t = compute_t(i, l, ct_x, ct_y);

	ct_z = HE_Context->hom_mult_opt(ct_z, ct_t);  // ct_z *= ct_t;
	*ret = HE_Context->hom_add(*ret, ct_z);		// ret += ct_z;

	return *ret;
}

/*************************************************************************************/
Mat_ZZ& compute_s(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y)
{
	assert(ct_x.size() > 0 && ct_y.size() > 0);
	if (j == 1)
	{
		Mat_ZZ *ret = new Mat_ZZ(ct_x[i]);
		*ret = HE_Context->hom_mult_opt(*ret, ct_y[i]);		// ret *= ct_y[i];
		*ret = HE_Context->hom_add(*ret, ct_y[i]);			// ret += ct_y[i];
		*ret = HE_Context->hom_add(*ret, (*ctxt_of_1));		// ret += *Mat_ZZ_of_1;
		return *ret;
	}

	int l;
	l = (j % 2 == 0) ? j / 2 : j / 2 + 1;

	Mat_ZZ *ret = new Mat_ZZ( compute_t(i + l, j - l, ct_x, ct_y) );
	Mat_ZZ ct_z = compute_z(i + l, j - l, ct_x, ct_y);
	Mat_ZZ ct_s = compute_s(i, l, ct_x, ct_y);

	ct_z = HE_Context->hom_mult_opt(ct_z, ct_s);	// ct_z *= ct_s;
	*ret = HE_Context->hom_add( *ret, ct_z);		    // ret += ct_z;

	return *ret;
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
		vt2[i] = HE_Context->hom_add(vt2[i], *ctxt_of_1);    // vt2[i] += ENC(1);
		vt2[i] = HE_Context->hom_mult_opt(vt2[i], b[i]);     // vt2[i] *= b[i];
		vt1[i] = HE_Context->hom_add(vt1[i], vt2[i]);        // vt1[i] += vt2[i];
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
