/*
 * Patch.cpp
 *
 *  Created on: 23 Dec, 2014
 *      Author: eeuser
 */

#include "Patch.h"

Patch::Patch() {
	pos3d = cv::Mat::zeros(4, 1, CV_64F);
	normal = cv::Mat::zeros(4, 1, CV_64F);
	consistentTex.clear();
	nonConsistentTex.clear();

	reProjectionsReady = false;
	reProjectionsColorsReady = false;

	spacialNNReady = false;
	optimalB=optimalG=optimalR=0;

}

Patch::~Patch() {
	// TODO Auto-generated destructor stub
}

const cv::Mat& Patch::getNormal() const {
	return normal;
}

void Patch::setNormal(double x, double y, double z, double w) {
	normal.at<double>(0) = x;
	normal.at<double>(1) = y;
	normal.at<double>(2) = z;
	normal.at<double>(3) = w;
}

const cv::Mat& Patch::getPos3d() const {
	return pos3d;
}

void Patch::setPos3d(double x, double y, double z, double w) {
	pos3d.at<double>(0) = x;
	pos3d.at<double>(1) = y;
	pos3d.at<double>(2) = z;
	pos3d.at<double>(3) = w;
}

float Patch::getPhotometricConsistencyScore() const {
	return photometricConsistencyScore;
}

void Patch::setPhotometricConsistencyScore(double photometricConsistencyScore) {
	this->photometricConsistencyScore = photometricConsistencyScore;
}

bool Patch::isReProjectionsReady() const {
	return reProjectionsReady;
}



void Patch::computeReprojections(vector<cv::Mat> ProjMat)
{
	if( reProjectionsReady )
	{
		cout<< "[Patch] Reprojections are already computed...NOT computing again\n";
		return;
	}

	consistentTexCords.reserve( consistentTex.size()  );
	for( int c=0; c<consistentTex.size() ; c++ )
	{
		int imIndx = consistentTex[c];
		cv::Mat rPt = ProjMat[imIndx] * pos3d;

		//convert homogeneous cord to normal cord and store it
		double tmpX = rPt.at<double>(0);
		double tmpY = rPt.at<double>(1);
		double tmpZ = rPt.at<double>(2);
		//cout << cv::Point2f(tmpX/tmpZ, tmpY/tmpZ)  << endl;
		consistentTexCords.push_back( cv::Point2f(tmpX/tmpZ, tmpY/tmpZ ) );
	}

	nonConsistentTexCords.reserve(nonConsistentTex.size());
	for( int n=0; n<nonConsistentTex.size() ; n++ )
	{
		int imIndx = nonConsistentTex[n];
		cv::Mat rPt = ProjMat[imIndx] * pos3d;

		//convert homogeneous cord to normal cord and store it
		double tmpX = rPt.at<double>(0);
		double tmpY = rPt.at<double>(1);
		double tmpZ = rPt.at<double>(2);
		//cout << cv::Point2f(tmpX/tmpZ, tmpY/tmpZ)  << endl;
		nonConsistentTexCords.push_back( cv::Point2f(tmpX/tmpZ, tmpY/tmpZ ) );

	}

	reProjectionsReady = true;
}

bool Patch::isReProjectionsColorsReady() const {
	return reProjectionsColorsReady;
}


// from the available reprojection co-ordinates get their colors
void Patch::getColorsAtReprojections( vector<cv::Mat> imgs )
{
	if( !reProjectionsReady )
	{
		cerr<< "[Patch] Reprojection co-ordinates not computed. Cannot compute their colors\nCall computeReprojections() first\n";
		return;
	}

	if( reProjectionsColorsReady )
	{
		cerr<< "[Patch] Colors already computed...NOT computing again\n";
		return;
	}

	avgB = avgG = avgR = 0; //also store the average color of each patch
	avgL = avga = avgb = 0; //averages of L a b
	consistentTexColors.reserve( consistentTex.size() );
	for( int i=0 ; i<consistentTex.size() ; i++ )
	{
		int imIndx = consistentTex[i];
		cv::Point2f cord = consistentTexCords[i];
		cv::Vec3b color = imgs[imIndx].at<cv::Vec3b>( (int)floor(cord.y), (int)floor(cord.x) );


		float tmpl, tmpa, tmpb;
		rgb2lab((float)color[2],(float)color[1],(float)color[0],  tmpl,tmpa,tmpb);
		avgL += tmpl;
		avga += tmpa;
		avgb += tmpb;

		avgB += (int)color[0];
		avgG += (int)color[1];
		avgR += (int)color[2];

		consistentTexColors.push_back( color );
	}

	//for( int i=0 ; i<consistentTexColors.size() ; i++ )
	//	cout<< consistentTexColors[i] << endl;

	//int tmo;
	//cin >> tmo;

	int len = consistentTex.size();


	//cout << avgB << " " << avgG << " " << avgR << endl;


	nonConsistentTexColors.reserve( nonConsistentTex.size() );
	for( int i=0 ; i<nonConsistentTex.size() ; i++ )
	{
		int imIndx = nonConsistentTex[i];
		cv::Point2f cord = nonConsistentTexCords[i];
		cv::Vec3b color = imgs[imIndx].at<cv::Vec3b>( (int)floor(cord.y), (int)floor(cord.x) );

		avgB += (int)color[0];
		avgG += (int)color[1];
		avgR += (int)color[2];
		nonConsistentTexColors.push_back( color );
	}

	int len2 = nonConsistentTex.size();
	avgB /= (len+len2);
	avgG /= (len+len2);
	avgR /= (len+len2);

	avgL /= len;
	avga /= len;
	avgb /= len;
	lab2rgb(avgL,avga,avgb,  xavgR,xavgG,xavgB);

	reProjectionsColorsReady = true;

}

void Patch::insertScribbles(vector<cv::Mat> gray3Channel)
{
	//cv::Vec3b col = cv::Vec3b( (uchar)avgB, (uchar)avgG, (uchar)avgR ); //simple set new colors as average
	cv::Vec3b col = cv::Vec3b( (uchar)optimalB, (uchar)optimalG, (uchar)optimalR );
	//cv::Vec3b col = cv::Vec3b( (uchar)xavgB, (uchar)xavgG, (uchar)xavgR );

	//cv::Vec3b col = cv::Vec3b( (uchar) this->getColor3dB(),(uchar) this->getColor3dG(),(uchar) this->getColor3dR() );

	//cv::Vec3b col = consistentTexColors[0];
	//cv::Vec3b col = cv::Vec3b( (uchar)0, (uchar)0, (uchar)0 );

	//cout<< "[Patch::insertScribbles] Optimal colors : "<< col << endl;


	for( int i=0 ; i<consistentTex.size() ; i++ )
	{
		int imIndx = consistentTex[i];
		gray3Channel[imIndx].at<cv::Vec3b>( consistentTexCords[i] ) = col;
	}

	for( int i=0 ; i<nonConsistentTex.size() ; i++ )
	{
		int imIndx = nonConsistentTex[i];
		//cout<< nonConsistentTexCords[i] << endl;
		if( nonConsistentTexCords[i].x >= 0  && nonConsistentTexCords[i].y >= 0 )
			gray3Channel[imIndx].at<cv::Vec3b>( nonConsistentTexCords[i] ) = col;
	}

}

bool Patch::isSpacialNnReady() const {
	return spacialNNReady;
}

void Patch::setNearestNeighbours(int n1, int n2, int n3, int n4)
{
	spacialNearestNeighbourIdx.reserve(4);
	spacialNearestNeighbourIdx.push_back(n1);
	spacialNearestNeighbourIdx.push_back(n2);
	spacialNearestNeighbourIdx.push_back(n3);
	spacialNearestNeighbourIdx.push_back(n4);

	spacialNNReady = true;
}

void Patch::setNearestNeighbours(ANNidxArray ary, ANNdistArray dist, int count)
{
	spacialNearestNeighbourIdx.reserve(count);
	spacialNearestNeighbourDist.reserve(count);
	for( int i=0 ; i<count ; i++ )
	{
		if( dist[i] > 0 ) { // only retain for +ve distances (avoid 0-distances)
			spacialNearestNeighbourIdx.push_back(ary[i]);
			spacialNearestNeighbourDist.push_back(dist[i]);
		}
	}

	spacialNNReady = true;
}



// using http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void Patch::rgb2lab( float R, float G, float B, float & l_s, float &a_s, float &b_s )
{
	float var_R = R/255.0;
	float var_G = G/255.0;
	float var_B = B/255.0;


	if ( var_R > 0.04045 ) var_R = pow( (( var_R + 0.055 ) / 1.055 ), 2.4 );
	else                   var_R = var_R / 12.92;
	if ( var_G > 0.04045 ) var_G = pow( ( ( var_G + 0.055 ) / 1.055 ), 2.4);
	else                   var_G = var_G / 12.92;
	if ( var_B > 0.04045 ) var_B = pow( ( ( var_B + 0.055 ) / 1.055 ), 2.4);
	else                   var_B = var_B / 12.92;

	var_R = var_R * 100.;
	var_G = var_G * 100.;
	var_B = var_B * 100.;

	//Observer. = 2°, Illuminant = D65
	float X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
	float Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
	float Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;


	float var_X = X / 95.047 ;         //ref_X =  95.047   Observer= 2°, Illuminant= D65
	float var_Y = Y / 100.000;          //ref_Y = 100.000
	float var_Z = Z / 108.883;          //ref_Z = 108.883

	if ( var_X > 0.008856 ) var_X = pow(var_X , ( 1./3. ) );
	else                    var_X = ( 7.787 * var_X ) + ( 16. / 116. );
	if ( var_Y > 0.008856 ) var_Y = pow(var_Y , ( 1./3. ));
	else                    var_Y = ( 7.787 * var_Y ) + ( 16. / 116. );
	if ( var_Z > 0.008856 ) var_Z = pow(var_Z , ( 1./3. ));
	else                    var_Z = ( 7.787 * var_Z ) + ( 16. / 116. );

	l_s = ( 116. * var_Y ) - 16.;
	a_s = 500. * ( var_X - var_Y );
	b_s = 200. * ( var_Y - var_Z );


}

//http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void Patch::lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B )
{
	float var_Y = ( l_s + 16. ) / 116.;
	float var_X = a_s / 500. + var_Y;
	float var_Z = var_Y - b_s / 200.;

	if ( pow(var_Y,3) > 0.008856 ) var_Y = pow(var_Y,3);
	else                      var_Y = ( var_Y - 16. / 116. ) / 7.787;
	if ( pow(var_X,3) > 0.008856 ) var_X = pow(var_X,3);
	else                      var_X = ( var_X - 16. / 116. ) / 7.787;
	if ( pow(var_Z,3) > 0.008856 ) var_Z = pow(var_Z,3);
	else                      var_Z = ( var_Z - 16. / 116. ) / 7.787;

	float X = 95.047 * var_X ;    //ref_X =  95.047     Observer= 2°, Illuminant= D65
	float Y = 100.000 * var_Y  ;   //ref_Y = 100.000
	float Z = 108.883 * var_Z ;    //ref_Z = 108.883


	var_X = X / 100. ;       //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
	var_Y = Y / 100. ;       //Y from 0 to 100.000
	var_Z = Z / 100. ;      //Z from 0 to 108.883

	float var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
	float var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
	float var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

	if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R , ( 1 / 2.4 ))  - 0.055;
	else                     var_R = 12.92 * var_R;
	if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G , ( 1 / 2.4 ) )  - 0.055;
	else                     var_G = 12.92 * var_G;
	if ( var_B > 0.0031308 ) var_B = 1.055 * pow( var_B , ( 1 / 2.4 ) ) - 0.055;
	else                     var_B = 12.92 * var_B;

	R = var_R * 255.;
	G = var_G * 255.;
	B = var_B * 255.;

}
