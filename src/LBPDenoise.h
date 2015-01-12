/*
 * LoopyBPDenoising.h
 *
 *  Created on: 27 Dec, 2014
 *      Author: eeuser
 */

#ifndef LOOPYBPDENOISING_H_
#define LOOPYBPDENOISING_H_

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <time.h>



//#define print(r) cout << r << endl;
#define print(r)


#define __OPTIMIZED_ONLY

/*
#ifdef __ORIGINAL_ONLY
#define opt(r)
#define optK(r)
#define org(r) r
#endif

#ifdef __OPTIMIZED_ONLY
#define opt(r) r ; x=y;
#define optK(r) r ;
#define org(r)
#endif
*/
/*
#define opt(r) r ; x=y;
//#define opt(r)


#define optK(r) r ;
//#define optK(r)

//#define org(r) r
#define org(r)
*/

using namespace std;


class Position
{
public:
	Position() { r=-1 ; c = -1; }
	Position( int r1, int c1 ) { r = r1; c = c1; }

	void set( int r1, int c1 ) { r = r1; c = c1; }
	int r, c;
};


class LoopyBPDenoising {
public:
	LoopyBPDenoising(string fileName );
	virtual ~LoopyBPDenoising();

	bool loadImage();
	const cv::Mat& getIm() const;

	bool denoise();

private:
	string imFileName;
	cv::Mat im;
	cv::Mat output;
	int nRows, nCols, nPixels;

	int nLabels;

	//incoming message-box for every pixel.
	// For every pixel it will be of size nLabel x 5 (nNeighbours)
	vector< cv::Mat > xl; //incoming messages from left
	vector< cv::Mat > xr; //      --"--       from right
	vector< cv::Mat > xt; //      --"--       from top
	vector< cv::Mat > xb; // 	 --"-- 		 from bottom
	vector< cv::Mat > belief; // the belief vector


	double param_tau, param_s, param_d;

	////// Felzenszwalb's Optimization Related ///////

	//kappa function
	// h_p(fp) = min_{f_p} (  h(f_p) + d  )
	double kappaRight( Position p );
	double passMessageRightOpt( Position p, Position q, int fq, double kappa );

	double kappaLeft( Position p );
	double passMessageLeftOpt( Position p, Position q, int fq, double kappa );

	double kappaTop( Position p );
	double passMessageTopOpt( Position p, Position q, int fq, double kappa );

	double kappaBottom( Position p );
	double passMessageBottomOpt( Position p, Position q, int fq, double kappa );
	/////////////////////////////////////////////////



	double dataCost( int y, int l ) { return std::min( (double)(abs(y-l)), param_tau) ; } // y is the intensity at this position
	//double dataCost( int y, int l ) { return std::min( (double)(abs(y-16*l)), param_tau) ; } // y is the intensity at this position

	//double smoothnessCost( int l, int l1 ) { int n=l-l1; return( -std::min((double)abs(n), param_d) );}
	double smoothnessCost( int l, int l1 ) { int n=l-l1; if( n==0) return 0.0; else return param_d;}

	void debug();

	void initMessages();

	void passMessageRight();
	void passMessageRight(Position p, Position q);
	double passMessageRight( Position p, Position q, int fq );


	void passMessageLeft();
	void passMessageLeft(Position p, Position q);
	double passMessageLeft( Position p, Position q, int fq );


	void passMessageTop();
	void passMessageTop(Position p, Position q);
	double passMessageTop( Position p, Position q, int fq );


	void passMessageBottom();
	void passMessageBottom(Position p, Position q);
	double passMessageBottom( Position p, Position q, int fq );

	void normalize( Position p, Position q, int fq, double factor, vector<cv::Mat> tU  );


	void computeBelief( );
	void computeBelief( Position q );
	void computeBelief( Position q, int fq);


	void makeOutputFromBelief();
	void makeOutputFromBelief( Position p );



	void writeMatToFile(cv::Mat& m, const char* filename);


	void print_elapsed_time( time_t start, time_t end) {
		double dif = difftime (end,start);
		printf ("Elasped time is %.2lf seconds.", dif );
	}


};

#endif /* LOOPYBPDENOISING_H_ */


// References :
// http://nghiaho.com/?page_id=1366
// file:///home/eeuser/Downloads/BeliefPropagation.pdf
// https://classes.soe.ucsc.edu/cmps290c/Spring04/proj/BPApp.pdf

// For speed up : http://www.cs.cornell.edu/~dph/papers/bp-cvpr.pdf
