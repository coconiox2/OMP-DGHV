#include "teste.h"

Flat_DGHV* HE_Context;
Mat_ZZ* ctxt_of_1;                  // encryption of 1 (constant)
int		t_bits;
ZZ *baza;


int main(int argc, char **argv)
{
    int lambda = 192;
    int int_size = 8;
    int N = 16;
    
    if(argc == 4)
    {
        lambda = atoi(argv[1]);
        int_size = atoi(argv[2]);
        N = atoi(argv[3]);
    }
    
    // test_max_no_encryptions();
    
    // cout << "\n\nSECURITATE COMPROMISA prin reducerea cheii publice la tau = lambda de la tau = O(gamma+lambda).";
    //  cout << endl << endl;
    
    // test_max_mult_depth();
    
    // test_schema_FDGHV();
    
    test_max_func(lambda, int_size, N);
    
    // cmp_HElib();
    
	// cout << " bitset pentru matricile GSW\n\n";
	
	// test_Flat_DGHV_simple();
    
    // test_FDGHV_max_depth();

	return 0;

}
