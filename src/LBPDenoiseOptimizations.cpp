/*
 * LBPDenoiseOptimizations.cpp
 *
 *  Created on: 30 Dec, 2014
 *      Author: eeuser
 */

#include "LBPDenoise.h"


/////////////////////// START OPTIMIZED FUNCTION (RIGHT) ///////////////////////////////
//kappa function
// h_p(fp) = min_{f_p} (  h(f_p) + d  )
double LoopyBPDenoising::kappaRight( Position p )
{
	//vector<double> tmpVec(nLabels);
	//tmpVec.clear();
	double miniVal=1.0e20;
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = dataCost(I_p,fp) ;
		thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		//tmpVec.push_back(thisC);
		if( miniVal > thisC )
			miniVal = thisC;
	}



	//double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	return (miniVal+param_d);
}

// for a given pixel, computes a component of the vector, optimized
double LoopyBPDenoising::passMessageRightOpt( Position p, Position q, int fq, double kappa )
{
	int I_q = (int)im.at<uchar>(p.r,p.c);
	double h_fq = dataCost(I_q,fq);
	h_fq += xt[fq].at<double>(p.r,p.c);
	h_fq += xb[fq].at<double>(p.r,p.c);
	h_fq += xl[fq].at<double>(p.r,p.c);
	h_fq += xr[fq].at<double>(p.r,p.c);


	double toSet = std::min( h_fq, kappa );
	xl[fq].at<double>(q.r,q.c) = toSet;
	return toSet;
}

/////////////////////// END OPTIMIZED FUNCTION (RIGHT)///////////////////////////////



/////////////////////// START OPTIMIZED FUNCTION (LEFT)///////////////////////////////
//kappa function
// h_p(fp) = min_{f_p} (  h(f_p) + d  )
double LoopyBPDenoising::kappaLeft( Position p )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = dataCost(I_p,fp);
		thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	return (minVal+param_d);
}

// for a given pixel, computes a component of the vector, optimized
double LoopyBPDenoising::passMessageLeftOpt( Position p, Position q, int fq, double kappa )
{
	int I_q = (int)im.at<uchar>(p.r,p.c);
	double h_fq = dataCost(I_q,fq);
	h_fq += xt[fq].at<double>(p.r,p.c);
	h_fq += xb[fq].at<double>(p.r,p.c);
	h_fq += xl[fq].at<double>(p.r,p.c);
	h_fq += xr[fq].at<double>(p.r,p.c);

	//cout<< h_fq << " "<< kappa << endl;
	double toSet = std::min( h_fq, kappa );
	xr[fq].at<double>(q.r,q.c) = toSet;
	return toSet;
}

/////////////////////// END OPTIMIZED FUNCTION (LEFT)///////////////////////////////



/////////////////////// START OPTIMIZED FUNCTION (TOP)///////////////////////////////
//kappa function
// h_p(fp) = min_{f_p} (  h(f_p) + d  )
double LoopyBPDenoising::kappaTop( Position p )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = dataCost(I_p,fp);
		thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	return (minVal+param_d);
}

// for a given pixel, computes a component of the vector, optimized
double LoopyBPDenoising::passMessageTopOpt( Position p, Position q, int fq, double kappa )
{
	int I_q = (int)im.at<uchar>(p.r,p.c);
	double h_fq = dataCost(I_q,fq);
	h_fq += xt[fq].at<double>(p.r,p.c);
	h_fq += xb[fq].at<double>(p.r,p.c);
	h_fq += xl[fq].at<double>(p.r,p.c);
	h_fq += xr[fq].at<double>(p.r,p.c);

	//cout<< h_fq << " "<< kappa << endl;
	double toSet = std::min( h_fq, kappa );
	xb[fq].at<double>(q.r,q.c) = toSet;
	return toSet;
}

/////////////////////// END OPTIMIZED FUNCTION (TOP)///////////////////////////////



/////////////////////// START OPTIMIZED FUNCTION (BOTTOM)///////////////////////////////
//kappa function
// h_p(fp) = min_{f_p} (  h(f_p) + d  )
double LoopyBPDenoising::kappaBottom( Position p )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = dataCost(I_p,fp);
		thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	return (minVal+param_d);
}

// for a given pixel, computes a component of the vector, optimized
double LoopyBPDenoising::passMessageBottomOpt( Position p, Position q, int fq, double kappa )
{
	int I_q = (int)im.at<uchar>(p.r,p.c);
	double h_fq = dataCost(I_q,fq);
	h_fq += xt[fq].at<double>(p.r,p.c);
	h_fq += xb[fq].at<double>(p.r,p.c);
	h_fq += xl[fq].at<double>(p.r,p.c);
	h_fq += xr[fq].at<double>(p.r,p.c);

	//cout<< h_fq << " "<< kappa << endl;
	double toSet = std::min( h_fq, kappa );
	xt[fq].at<double>(q.r,q.c) = toSet;
	return toSet;
}

/////////////////////// END OPTIMIZED FUNCTION (BOTTOM)///////////////////////////////

