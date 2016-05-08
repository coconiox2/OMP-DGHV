#pragma once

#include <iostream>
#include <vector>
#include "../Code/Flat_DGHV.h"

using namespace std;

/* integer bitwise HE */
class HE_Integer
{
	int	m_t_bits;	// working length
	const Flat_DGHV& he_context;
public:
	HE_Integer(const Flat_DGHV& he_context, int t_bits);

	void	encryptBit(Mat_ZZ& ct, int bit);
	void	encryptIntValue(vector<Mat_ZZ>& vct, int val);

	int		decryptBit(const Mat_ZZ &ct);
	int 	decryptIntValue(const vector<Mat_ZZ>& vct);

	void	encryptIntVector(vector<vector<Mat_ZZ> >& vvct, const vector<int>& values);
};

