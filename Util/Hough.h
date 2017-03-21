#ifndef Hough_H
#define Hough_H
#include "const.h"

void HoughTransform( std::vector<Double_t> &f_X,   std::vector<Double_t> &f_Y,
		     std::vector<Double_t> &f_afX, std::vector<Double_t> &f_afY );
void HoughFit_One( TH2D* f_hist, Double_t& f_par0, Double_t& f_par1 );

void GetPol1Residual (TF1 *f_func, std::vector<Double_t> &f_X, std::vector<Double_t> &f_Y, std::vector<Double_t>& f_resi );
void GetPol1ResidualY(TF1 *f_func, std::vector<Double_t> &f_X, std::vector<Double_t> &f_Y, std::vector<Double_t>& f_resi );

#endif
