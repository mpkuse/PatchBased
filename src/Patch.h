/*
 * Patch.h
 *
 *  Created on: 23 Dec, 2014
 *      Author: eeuser
 */

#ifndef PATCH_H_
#define PATCH_H_

#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>


#include <ANN/ANN.h>

using namespace std;

class Patch {
public:
	Patch();
	virtual ~Patch();
	const cv::Mat& getNormal() const;
	void setNormal(double x, double y, double z, double w);
	const cv::Mat& getPos3d() const;
	void setPos3d(double x, double y, double z, double w);
	float getPhotometricConsistencyScore() const;
	void setPhotometricConsistencyScore(double photometricConsistencyScore);

	void setDebug(float d1, float d2) {
		this->debug1 = d1;
		this->debug2 = d2;
	}

	bool isReProjectionsReady() const;

	bool isReProjectionsColorsReady() const;

	vector<int> consistentTex; //image index list(0 based) in which the point is visible and texture agrees also
	vector<int> nonConsistentTex; //image index list in which the point is visible but texture may not agree

	vector<cv::Point2f> consistentTexCords; //re-projected cord position in each visible view.
	vector<cv::Point2f> nonConsistentTexCords;

	vector<cv::Vec3b> consistentTexColors; //colors at corresponding consistentTexCords
	vector<cv::Vec3b> nonConsistentTexColors;

	// given all the projection matrices computes the reprojection of this patch on all visible views
	void computeReprojections( vector<cv::Mat> ProjMat);

	// from the available reprojection co-ordinates get their colors
	void getColorsAtReprojections( vector<cv::Mat> imgs );

	int avgR, avgG, avgB;
	float avgL, avga, avgb;
	float xavgR, xavgG, xavgB;

	void insertScribbles( vector<cv::Mat> gray3Channel );
	bool isSpacialNnReady() const;

	// to compute nn with ANN and store the index here. Interested in 4 nearest neighbours.
	//    memory for this is reserved in the constructor.
	vector<int> spacialNearestNeighbourIdx;
	vector<double> spacialNearestNeighbourDist; //corresponding euclidean distances
	void setNearestNeighbours( int n1, int n2, int n3, int n4 ); // pass as arg the indices of nearest-neighs
	void setNearestNeighbours( ANNidxArray ary, ANNdistArray dist, int count);

	int getOptimalB() const {
		return optimalB;
	}

	void setOptimalB(int optimalB) {
		this->optimalB = optimalB;
	}

	int getOptimalG() const {
		return optimalG;
	}

	void setOptimalG(int optimalG) {
		this->optimalG = optimalG;
	}

	int getOptimalR() const {
		return optimalR;
	}

	void setOptimalR(int optimalR) {
		this->optimalR = optimalR;
	}

	int getColor3dB() const {
		return color3dB;
	}

	void setColor3dB(int color3dB) {
		this->color3dB = color3dB;
	}

	int getColor3dG() const {
		return color3dG;
	}

	void setColor3dG(int color3dG) {
		this->color3dG = color3dG;
	}

	int getColor3dR() const {
		return color3dR;
	}

	void setColor3dR(int color3dR) {
		this->color3dR = color3dR;
	}

private:
	cv::Mat pos3d; //3d position in homogeneous cords
	cv::Mat normal; //normal in homogenerous cords. note normal.w will be zero (since it denotes a direction)
	double photometricConsistencyScore;
	float debug1, debug2;

	bool reProjectionsReady;
	bool reProjectionsColorsReady;

	bool spacialNNReady;


	int optimalR, optimalG, optimalB;


	int color3dR, color3dG, color3dB;

	// using http://www.easyrgb.com/index.php?X=MATH&H=01#text1
	void rgb2lab( float R, float G, float B, float & l_s, float &a_s, float &b_s );


	//http://www.easyrgb.com/index.php?X=MATH&H=01#text1
	void lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B );

};

#endif /* PATCH_H_ */
