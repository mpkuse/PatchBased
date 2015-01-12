/*
 * ApproxNN.cpp
 *
 *  Created on: 2 Jan, 2015
 *      Author: eeuser
 */

#include "Patch.h"
#include "Patches.h"
/********************
 *  TEST code for ANN
 ***********************
void printPt(ANNpoint p)			// print point
{
	cout << "(" << p[0];
	for (int i = 1; i < DIMS; i++) {
		cout << ", " << p[i];
	}
	cout << ")\n";
}

// reads the datapoint from `class Patch` and stores into ANNpoint.
void readDataPt( ANNpoint pt, const Patch & patch )
{
	cv::Mat matPt = patch.getPos3d();
	for( int i=0 ; i<DIMS ; i++ )
		pt[i] = matPt.at<double>(i);
}


void testANN( Patches & p )
{
	ANNpointArray dataPts = annAllocPts(p.getPatchCount(), DIMS);
	ANNpoint queryPt = annAllocPt(DIMS);
	//queryPt[0] = -1.08199;
	//queryPt[1] = 0.894556;
	//queryPt[2] = -3.61927;
	readDataPt(queryPt, p.getPatchAt(0));
	printPt(queryPt);


	ANNidxArray nnIdx = new ANNidx[nNN];
	ANNdistArray dists = new ANNdist[nNN];

	cout<< "[testANN] Fill ANNpointArray\n";
	for( int i=0 ; i<p.getPatchCount() ; i++ )
	{
		readDataPt(dataPts[i], p.getPatchAt(i));
		//printPt(dataPts[i]);
	}

	cout<< "[testAnn] Building KD-tree\n";
	ANNkd_tree * kdTree = new ANNkd_tree(  dataPts, p.getPatchCount(), DIMS );


	cout<< "[testAnn] Doing annkSearch\n";
	kdTree->annkSearch(queryPt,nNN, nnIdx, dists, 0.0 );

	for( int i=1 ; i<nNN ; i++ )
	{
		cout<< "Nearest Nei index=" << nnIdx[i] << " ;;; dist= "<< sqrt(dists[i]) << endl;
	}


	//ofstream ff("file.txt");
	//kdTree->Dump(ANNtrue, ff);
	//ff.close();

	delete [] nnIdx;
	delete [] dists;
	annClose();
}

*********************END OF TEST CODE*************/
////////////////////////////////////////////////////////////////////////////////////////
///////////////////			CLASS RELATED    	///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////



// nearest neighbour graph related
void Patches::computeSpacialNearestNeighbourGraph()
{
	cout<< "[Patches/ANN] computeSpacialNearestNeighbourGraph:\n";

	// create ANNpointArray ;;; DIMS = 3
	ANNpointArray dataPts = annAllocPts(patchVector.size(), DIMS);
	ANNpoint queryPt = annAllocPt(DIMS);

	ANNidxArray nnIdx = new ANNidx[nNN];
	ANNdistArray dists = new ANNdist[nNN];


	// load data into ANNpointArray
	for( int i=0 ; i<patchVector.size() ; i++ )
	{
		readDataPt(dataPts[i], patchVector[i]);
		//printPt(dataPts[i]);
	}

	// create KD-tree;;; nNN = 5
	cout<< "[Patches/ANN] Building KD-tree\n";
	ANNkd_tree * kdTree = new ANNkd_tree(  dataPts, patchVector.size(), DIMS );

	// for each query, find the nearest neighbours.
	cout<< "[Patches/ANN] Finding "<< nNN << " approximate nearest neighbour of every Patch\n";
	for( int q=0 ; q<patchVector.size() ; q++ )
	{
		readDataPt(queryPt, patchVector[q]);
		kdTree->annkSearch(queryPt,nNN, nnIdx, dists, 0.0 );

		// note not putting nnIdx since it will correspond to the queryPt. Each pt
		// is nearest to itself.
		//patchVector[q].setNearestNeighbours(nnIdx[1], nnIdx[2], nnIdx[3], nnIdx[4]);
		patchVector[q].setNearestNeighbours( nnIdx, dists, nNN );
	}


//	cv::Vec3b vf = patchVector[0].consistentTexColors[0];
//	cout<< "xxxx"<< (int)vf[0] << endl;
//
//	cout<< patchVector[1].spacialNearestNeighbourIdx.size() << endl;
//	cout<< patchVector[1].spacialNearestNeighbourIdx[0] << endl;
//	cout<< patchVector[1].spacialNearestNeighbourIdx[1] << endl;
//	cout<< patchVector[1].spacialNearestNeighbourIdx[2] << endl;
//	cout<< patchVector[1].spacialNearestNeighbourIdx[3] << endl;

	// deallocate memory
	delete [] nnIdx;
	delete [] dists;
	annClose();
}


void Patches::printPt(ANNpoint p)			// print point
{
	cout << "(" << p[0];
	for (int i = 1; i < DIMS; i++) {
		cout << ", " << p[i];
	}
	cout << ")\n";
}

// get jth nearest neighbour of patchVector[i];
const Patch& Patches::getNearestNeighbourPatchOf(int i, int j)
{
	if( patchVector[i].isSpacialNnReady() == false )
	{
		cerr<< "[Patches] Nearest neighbours are not computed\n...Quiting\n";
		exit(1);
	}

	if( j >= patchVector[i].spacialNearestNeighbourIdx.size() )
	{
		cerr<< "[Patches] j is larger than spacialNearestNeighbourIdx.size()\n";
		exit(1);
	}

	int returnIndx = patchVector[i].spacialNearestNeighbourIdx[j];
	return patchVector[ returnIndx ];

}

// reads the datapoint from `class Patch` and stores into ANNpoint.
void Patches::readDataPt( ANNpoint pt, const Patch & patch )
{
	cv::Mat matPt = patch.getPos3d();
	for( int i=0 ; i<DIMS ; i++ )
		pt[i] = matPt.at<double>(i);
}
