#include "HE_Integer.h"

/*************************************************************************************/
HE_Integer::HE_Integer(const Flat_DGHV& context, int t_bits)
	: he_context(context), m_t_bits(t_bits) {}

/*************************************************************************************/
void HE_Integer::encryptBit(Mat_ZZ& ct, int bit)
{
	ct = he_context.encrypt(bit);
}

void HE_Integer::encryptIntValue(vector<Mat_ZZ>& vct, int val)
{
	vct.clear();
	int bit = 0;
	do
	{
		bit = val % 2;
		Mat_ZZ ct = he_context.encrypt(bit);
		vct.push_back(ct);
		val /= 2;
	} while (val != 0);

	for (int i = vct.size(); i < m_t_bits; i++)
	{
		Mat_ZZ ct = he_context.encrypt(0);
		vct.push_back(ct);
	}
}

int	HE_Integer::decryptBit(const Mat_ZZ &ct)
{
	return he_context.decrypt(ct);
}

int HE_Integer::decryptIntValue(const vector<Mat_ZZ>& vct)
{
	int value = 0;
	double power_of_two = 1;

	for (int i = 0; i < vct.size(); i++)
	{
		value += power_of_two * he_context.decrypt(vct[i]);
		power_of_two *= 2;
	}

	return value;
}

void HE_Integer::encryptIntVector(vector<vector<Mat_ZZ> >& vvct, const vector<int>& values)
{
	vvct.clear();
	vvct.reserve(values.size());
	for (int i = 0; i < values.size(); i++)
	{
		vector<Mat_ZZ> encrypted_integer;
		encryptIntValue(encrypted_integer, values[i]);
		vvct.push_back(encrypted_integer);
	}
}


