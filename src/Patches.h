/*
 * Patches.h
 *
 *  Created on: 23 Dec, 2014
 *      Author: eeuser
 */

#ifndef PATCHES_H_
#define PATCHES_H_


#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;


#include "Patch.h"
#include <ANN/ANN.h>

#include "Potts.h"


class Patches {
public:
	Patches();
	virtual ~Patches();
	bool readPatchFile( );
	bool readImages( );
	bool showAllImages( );
	bool readProjectionMatrices( );
	const string& getPMVSDir() const;
	void setPMVSDir(const string& pmvsDir);

	const vector<cv::Mat>& getImgs() const;
	const Patch& getPatchAt(int indx) const;
	const vector<cv::Mat>& getProjMat() const;

	int getPatchCount() { return patchVector.size(); }

	void computeAllReprojections();
	void computeAllReprojectionsColors();

	void readPatchColors();

	void showAllImagesWithReprojections();

	void write4Colorization(); //write all the images in gray, images with scribbles


	// nearest neighbour graph related
	void computeSpacialNearestNeighbourGraph();


#define DIMS 3
#define nNN 5 //should be greater than or equal to 5
	// get jth nearest neighbour of patchVector[i];
	const Patch & getNearestNeighbourPatchOf( int i, int j );


	/////////////////////////////////////////////////////////
	//////////// MAX FLOW OPTIMIZATION RELATED //////////////
	/////////////////////////////////////////////////////////
	bool synchronizeColors();


	static void rgb2lab( float r, float g, float b, float & l_s, float &a_s, float & b_s );
	static void lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B );

	//////////////// END MAX FLOW RELATED ///////////////////

private:
	string pmvsDir;
	int nPatchs;
	vector<Patch> patchVector; // vector of patchs

	vector<cv::Mat> imgs; //all images
	vector<cv::Mat> projMat; //projection matrices

	void showThisPatch( Patch ch );

	void printPt(ANNpoint p);			// print point
	void readDataPt( ANNpoint pt, const Patch & patch );





	/////////////////////////////////////////////////////////
	//////////// MAX FLOW OPTIMIZATION RELATED //////////////
	/////////////////////////////////////////////////////////

	void addAllEdges(Potts * p);
	int getEdgeLambda( cv::Vec3b nodeAColor, cv::Vec3b nodeBColor );


	// given node i and label k
	int setupDataCost( Potts * p, int nLabels, int code );

	// computes the datacost of `patch` wrt to label k
	int dCostB( Patch P, int k );
	int dCostG( Patch P, int k );
	int dCostR( Patch P, int k );


	void print_elapsed_time( time_t start, time_t end) {
		double dif = difftime (end,start);
		printf ("Elasped time is %.2lf seconds\n", dif );
	}

	void collectOptimalSolution(Potts * p, int code);

	#define CODE_BLUE_CHANNEL (1)
	#define CODE_GREEN_CHANNEL (2)
	#define CODE_RED_CHANNEL (3)

	//////////////// END MAX FLOW RELATED ///////////////////

};

#endif /* PATCHES_H_ */
