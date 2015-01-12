/*
 * LoopyBPDenoising.cpp
 *
 *  Created on: 27 Dec, 2014
 *      Author: eeuser
 */

#include "LBPDenoise.h"

LoopyBPDenoising::LoopyBPDenoising( string fName ) {

	this->imFileName = fName;

	nLabels = 256;
	//nLabels = 16;


	param_tau = 100.0,
	param_s = 1.0;
	param_d = 20.0;

	loadImage();
}

LoopyBPDenoising::~LoopyBPDenoising() {
	// TODO Auto-generated destructor stub
}

bool LoopyBPDenoising::loadImage()
{
	if( !im.empty() )
	{
		cout<< "[LoopyBPDenoising] Image already loaded...Not loading again\n";
		return true;
	}

	im = cv::imread( imFileName.c_str(), CV_LOAD_IMAGE_GRAYSCALE );

	if( im.empty() )
	{
		cerr<< "[LoopyBPDenoising] Cannot load image : "<< imFileName.c_str() << endl;
		return false;
	}

	nRows = im.rows;
	nCols = im.cols;
	nPixels = nRows * nCols;
	return true;
}

const cv::Mat& LoopyBPDenoising::getIm() const {
	return im;
}

bool LoopyBPDenoising::denoise()
{
	if( im.empty() )
	{
		cerr<< "[LoopyBPDenoising] No image loaded. Load image before calling this function\n";
		return false;
	}


	initMessages();


	time_t start, end;

	for( int itr=0 ; itr< 20 ; itr++ )
	{

		debug();

		cout<< "Iteration "<< itr << endl;
		time(&start);
		passMessageRight(); cout<< "\nDone passing on the right\n";
		time(&end);
		print_elapsed_time(start,end);

		time(&start);
		passMessageBottom();  cout<< "\nDone passing on the down\n";
		time(&end);
		print_elapsed_time(start,end);

		time(&start);
		passMessageLeft();  cout<< "\nDone passing on the left\n";
		time(&end);
		print_elapsed_time(start,end);

		time(&start);
		passMessageTop();    cout<< "\nDone passing on the up\n";
		time(&end);
		print_elapsed_time(start,end);




		computeBelief(); cout<< "Done computing belief\n";
		makeOutputFromBelief(); cout<< "create output labels from belief\n";
		char oName[20]; sprintf( oName, "out%02d.png", itr );
		cv::imwrite(oName, output ); cout<< "write output to PNG file\n";
	}






	for( int i=0 ; i<nLabels ; i++ )
	{
		char fName[50];
		sprintf( fName, "msg/msg_%03d.txt", i );
		writeMatToFile( belief[i], fName );
	}


	return true;

}

void LoopyBPDenoising::debug()
{
	for( int fp=0 ; fp<nLabels ; fp++ )
	cout << setw(6)<< xr[fp].at<double>(20,30) << " ";
	cout<< "\n";
	for( int fp=0 ; fp<nLabels ; fp++ )
		cout <<setw(6)<< xl[fp].at<double>(20,30) << " ";
	cout<< "\n";
	for( int fp=0 ; fp<nLabels ; fp++ )
		cout <<setw(6)<< xt[fp].at<double>(20,30) << " ";
	cout<< "\n";
	for( int fp=0 ; fp<nLabels ; fp++ )
		cout <<setw(6)<< xb[fp].at<double>(20,30) << " ";
	cout<< "\n";
}

void LoopyBPDenoising::initMessages()
{
	cv::Mat z = cv::Mat::zeros( nRows, nCols, CV_64F );
	//z = z + 1;

	// have z stacked up nLabelsx5 times
	xl.reserve(nLabels);
	xr.reserve(nLabels);
	xb.reserve(nLabels);
	xt.reserve(nLabels);
	belief.reserve(nLabels);
	for( int i=0 ; i<nLabels ; i++ )
	{
		xl.push_back( z.clone() );
		xr.push_back( z.clone() );
		xb.push_back( z.clone() );
		xt.push_back( z.clone() );
		belief.push_back( z.clone() );
	}
}


void LoopyBPDenoising::normalize( Position p, Position q, int fq, double factor, vector<cv::Mat> tU )
{
	tU[fq].at<double>(q.r,q.c) -= factor;
}



//////////////////// RIGHT ////////////////
// all pixels
void LoopyBPDenoising::passMessageRight()
{
	for( int r=1 ; r<(nRows-1) ; r++ )
	{cout<< "."<< flush;
		for( int c=1 ; c<(nCols-1) ; c++ )
		{
			Position p(r,c);
			Position q(r,c+1);
			passMessageRight(p,q);
		}
	}
}


// for a given pixel, computes the vector of messages
void LoopyBPDenoising::passMessageRight( Position p, Position q )
{
#ifdef __OPTIMIZED_ONLY
	double kappa = kappaRight(p);
#endif

	double minVal = 1.0e30;
	double x;
	for( int fq=0 ; fq<nLabels ; fq++ ){
#ifndef  __OPTIMIZED_ONLY
		x= passMessageRight(p,q, fq);
#endif

#ifdef __OPTIMIZED_ONLY
		x=passMessageRightOpt(p,q,fq,kappa);
#endif

		if( minVal > x )
			minVal=x;
	}

	// normalize
	double factor = minVal;
	for( int fq=0 ; fq<nLabels ; fq++ )
			normalize(p,q, fq, factor, xl) ;

}




// for a given pixel, computes a component of the vector
double LoopyBPDenoising::passMessageRight( Position p, Position q, int fq )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = smoothnessCost(fp,fq) + dataCost(I_p,fp);
		thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		//thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	xl[fq].at<double>(q.r,q.c) = minVal;
	return minVal;

}

//////////////////// LEFT ////////////////
// all pixels
void LoopyBPDenoising::passMessageLeft()
{
	for( int r=1 ; r<(nRows-1) ; r++ )
	{cout<< "."<< flush;
		for( int c=1 ; c<(nCols-1) ; c++ )
		{
			Position p(r,c);
			Position q(r,c-1);
			passMessageLeft(p,q);
		}
	}
}


// for a given pixel, computes the vector of messages
void LoopyBPDenoising::passMessageLeft( Position p, Position q )
{
#ifdef __OPTIMIZED_ONLY
	double kappa = kappaLeft(p);
#endif

	double minVal = 1.0e30;
	double x;
	for( int fq=0 ; fq<nLabels ; fq++ ){
#ifndef  __OPTIMIZED_ONLY
		x= passMessageLeft(p,q, fq);
#endif

#ifdef __OPTIMIZED_ONLY
		x= passMessageLeftOpt(p,q, fq, kappa);
#endif

		if( minVal > x )
			minVal=x;
	}

	// normalize
	double factor = minVal;
	for( int fq=0 ; fq<nLabels ; fq++ )
			normalize(p,q, fq, factor, xr) ;
}





// for a given pixel, computes a component of the vector
double LoopyBPDenoising::passMessageLeft( Position p, Position q, int fq )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = smoothnessCost(fp,fq) + dataCost(I_p,fp);
		thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		//thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	xr[fq].at<double>(q.r,q.c) = minVal;
	return minVal;

}


//////////////////// Top ////////////////
// all pixels
void LoopyBPDenoising::passMessageTop()
{
	for( int r=1 ; r<(nRows-1) ; r++ )
	{cout<< "."<< flush;
		for( int c=1 ; c<(nCols-1) ; c++ )
		{
			Position p(r,c);
			Position q(r-1,c);
			passMessageTop(p,q);
		}
	}
}


// for a given pixel, computes the vector of messages
void LoopyBPDenoising::passMessageTop( Position p, Position q )
{
#ifdef __OPTIMIZED_ONLY
	double kappa = kappaTop(p);
#endif

	double minVal = 1.0e30;
	double x;
	for( int fq=0 ; fq<nLabels ; fq++ ) {
#ifndef  __OPTIMIZED_ONLY
		x=passMessageTop(p,q, fq);
#endif

#ifdef __OPTIMIZED_ONLY
		x=passMessageTopOpt(p,q,fq,kappa);
#endif

		if( minVal > x )
			minVal=x;
	}

	// normalize
	double factor = minVal;
	for( int fq=0 ; fq<nLabels ; fq++ )
			normalize(p,q, fq, factor, xb) ;
}





// for a given pixel, computes a component of the vector
double LoopyBPDenoising::passMessageTop( Position p, Position q, int fq )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = smoothnessCost(fp,fq) + dataCost(I_p,fp);
		//thisC += xt[fp].at<double>(p.r,p.c);
		thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	xb[fq].at<double>(q.r,q.c) = minVal;
	return minVal;

}




//////////////////// Bottom ////////////////
// all pixels
void LoopyBPDenoising::passMessageBottom()
{
	for( int r=1 ; r<(nRows-1) ; r++ )
	{cout<< "."<< flush;
		for( int c=1 ; c<(nCols-1) ; c++ )
		{
			Position p(r,c);
			Position q(r+1,c);
			passMessageBottom(p,q);
		}
	}
}


// for a given pixel, computes the vector of messages
void LoopyBPDenoising::passMessageBottom( Position p, Position q )
{
#ifdef __OPTIMIZED_ONLY
	double kappa = kappaBottom(p);
#endif

	double minVal = 1.0e30;
	double x;
	for( int fq=0 ; fq<nLabels ; fq++ ) {
#ifndef  __OPTIMIZED_ONLY
		x=passMessageBottom(p,q, fq);
#endif

#ifdef __OPTIMIZED_ONLY
		x=passMessageBottomOpt(p,q,fq,kappa);
#endif

		if( minVal > x )
			minVal=x;
	}

	// normalize
	double factor = minVal;
	for( int fq=0 ; fq<nLabels ; fq++ )
			normalize(p,q, fq, factor, xt) ;
}




// for a given pixel, computes a component of the vector
double LoopyBPDenoising::passMessageBottom( Position p, Position q, int fq )
{
	vector<double> tmpVec(nLabels);
	tmpVec.clear();
	for( int fp=0 ; fp<nLabels ; fp++ )
	{
		int I_p = (int)im.at<uchar>(p.r,p.c);
		double thisC = smoothnessCost(fp,fq) + dataCost(I_p,fp);
		thisC += xt[fp].at<double>(p.r,p.c);
		//thisC += xb[fp].at<double>(p.r,p.c);
		thisC += xl[fp].at<double>(p.r,p.c);
		thisC += xr[fp].at<double>(p.r,p.c);

		tmpVec.push_back(thisC);
	}

	double minVal = *(min_element( tmpVec.begin(), tmpVec.end() ));

	xt[fq].at<double>(q.r,q.c) = minVal;
	return minVal;

}


//////////////////// COMPUTING BELIEF //////////////////////
// for all pixels
void LoopyBPDenoising::computeBelief( )
{
	for( int r=1 ; r<(nRows-1) ; r++ )
	{
		for( int c=1 ; c<(nCols-1) ; c++ )
		{
			Position q(r,c);
			computeBelief( q );
		}
	}
}


// compute the belief vector at position q
void LoopyBPDenoising::computeBelief( Position q )
{
	for( int fq=0 ; fq<nLabels ; fq++ )
		computeBelief( q, fq );
}


// compute belief at fq^{th} component given position q,
void LoopyBPDenoising::computeBelief( Position q, int fq)
{
	int I_q = (int)im.at<uchar>(q.r,q.c);
	 double thisBelief = dataCost( I_q, fq );
	 thisBelief += xl[fq].at<double>(q.r,q.c);
	 thisBelief += xr[fq].at<double>(q.r,q.c);
	 thisBelief += xt[fq].at<double>(q.r,q.c);
	 thisBelief += xb[fq].at<double>(q.r,q.c);


	 belief[fq].at<double>(q.r,q.c) = thisBelief;
}


void LoopyBPDenoising::makeOutputFromBelief( )
{
	output = cv::Mat::zeros(nRows, nCols, im.type() );

	for( int r=1 ; r<(nRows-1) ; r++ )
	{
		for( int c=1 ; c<(nCols-1) ; c++ )
		{
			Position p(r,c);
			makeOutputFromBelief( p );
		}
	}
}


void LoopyBPDenoising::makeOutputFromBelief( Position p )
{
	double minVal=1e100;
	int minIndx = -1;
	for( int l=0 ; l<nLabels ; l++ )
	{
		double b = belief[l].at<double>(p.r,p.c);
		if( minVal > b )
		{
			minVal = b;
			minIndx = l;
		}
	}
	output.at<uchar>(p.r,p.c) = (uchar)(minIndx);
	//output.at<uchar>(p.r,p.c) = (uchar)(minIndx * 16);
}


void LoopyBPDenoising::writeMatToFile(cv::Mat& m, const char* filename)
{
    ofstream fout(filename,ifstream::trunc);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }

    for(int i=0; i<m.rows; i++)
    {
        for(int j=0; j<m.cols; j++)
        {
            fout<< (double)m.at<double>(i,j)<<"\t";
        }
        fout<<endl;
    }

    fout.close();
}
