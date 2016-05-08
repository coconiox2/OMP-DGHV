#pragma once

#include <iostream>
#include <vector>
#include "../Code/Flat_DGHV.h"


/*************************************************************************************/
/* these are used to evaluate the "comparison" */
Mat_ZZ& compute_z(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y);
Mat_ZZ& compute_t(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y);
Mat_ZZ& compute_s(int i, int j, vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y);

/*************************************************************************************/
/* comparisons evaluations */
// void evaluate_X_gt_Y(vector<Mat_ZZ>& ct_x, vector<Mat_ZZ>& ct_y, int t_bits);
// void evaluate_X_ge_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits);
// void evaluate_X_eq_Y(vector<BigPoly>& ct_x, vector<BigPoly>& ct_y, int t_bits);

/*************************************************************************************/
/* these are used to evaluate the "maximum" */
vector<Mat_ZZ> select(Mat_ZZ& c, vector<Mat_ZZ>& a, vector<Mat_ZZ>& b);
vector<Mat_ZZ> getmax(vector<vector<Mat_ZZ> >& vvct, int start, int n); // A tree approach with a small nr. of levels consumed
int gmax(vector<int>& v, int start, int n);