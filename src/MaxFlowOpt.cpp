/*
 * MaxFlowOpt.cpp
 *
 *  Created on: 4 Jan, 2015
 *      Author: eeuser
 *
 *      Implements the functions related to computation of optimum color of each patch
 *
 *      These functions will be private functions of `class Patches`
 */


#include "Patches.h"
#include "Potts.h"


bool Patches::synchronizeColors()
{
	cout << "[Patches/Max-flow] Synchronize Patch Colors with max-flow\n";

	// [ Potts model, parametric maxflow and k-submodular functions ]
	// http://pub.ist.ac.at/~vnk/software.html

	//check for readyness of nearest-neighbours & colors & cords
	for( int i=0 ; i<patchVector.size() ; i++)
	{
		if( patchVector[i].isReProjectionsColorsReady()==false
									|| patchVector[i].isSpacialNnReady()==false)
		{
			cerr<< "[Patches/Max-flow] Either of reprojection colors and nearest-neighbours not ready.\n"
					"\tPlease call the computeSpacialNearestNeighbourGraph() before calling this\n";
			return false;
		}
	}

	const int nodeCount = patchVector.size();
	const int edgeCount = patchVector.size()*nNN;
	const int labelCount = 256;
	const int nIterations = 25;

	/// Create Solver
	Potts* p[3];
	p[0] = new Potts( nodeCount, edgeCount, labelCount, NULL );
	p[1] = new Potts( nodeCount, edgeCount, labelCount, NULL );
	p[2] = new Potts( nodeCount, edgeCount, labelCount, NULL );


	/// Add edges
	cout<< "[Patches/Max-flow] Adding graph edges from nearest-neighbour graph\n";
	addAllEdges(p[0]);
	addAllEdges(p[1]);
	addAllEdges(p[2]);

	// for each channel
	int channel[] = { CODE_BLUE_CHANNEL, CODE_GREEN_CHANNEL, CODE_RED_CHANNEL };

	for( int ch=0 ; ch<=2 ; ch++ )
	{
		cout<< "-----------\nProcessing Channel : "<< ch << endl;
		cout<< "-----------\n";

		// Set datacost at every node. Each node has costs at all labels
		cout<< "[Patches/Max-flow] Setting up datacost at every node\n";
		setupDataCost(p[ch], labelCount, channel[ch]);



		// Iterative Solver
		cout<< "[Patches/Max-flow] .....Solving.........\n";
		time_t start, end;
		time(&start);
		int c_star = (p[ch])->Solve( nIterations, false );
		time(&end);
		print_elapsed_time(start,end);
		cout<< "Optimal Value : "<< c_star << endl;


		// Collect solution
		collectOptimalSolution(p[ch], channel[ch]);

		cout<< "[Patches/Max-flow] Number of unlabeled nodes : "<< p[ch]->GetUnlabeledNum() << endl;

		delete p[ch];
	}

	return true;

}

int Patches::getEdgeLambda( cv::Vec3b A, cv::Vec3b B )
{
	float UA = -.1471*(float) A[2] - .28888*(float) A[1] + .436*(float) A[0];
	float UB = -.1471*(float) B[2] - .28888*(float) B[1] + .436*(float) B[0];

	float VA = .615*(float) A[2] - .5149*(float) A[1] - .1*(float) A[0];
	float VB = .615*(float) B[2] - .5149*(float) B[1] - .1*(float) B[0];


	float dist = (UA-UB)*(UA-UB) + (VA-VB)*(VA-VB);

	if( dist > 40 )// different colors
		return 5;
	else
		return 30;
}


// adds all the edges. It has to be already ensured that the nearest neighbour graph is ready
void Patches::addAllEdges(Potts * p)
{
	for( int pa=0 ; pa<patchVector.size() ; pa++ )
	{
		int numNn = patchVector[pa].spacialNearestNeighbourIdx.size();
		for( int n=0 ; n<numNn ; n++ )
		{
			int neihIndx = patchVector[pa].spacialNearestNeighbourIdx[n];
			int nodeA = std::min( pa, neihIndx );
			int nodeB = std::max( pa, neihIndx );
			//int lambda = 20;
			int lambda = getEdgeLambda( patchVector[nodeA].consistentTexColors[0], patchVector[nodeB].consistentTexColors[0] );

			p->AddEdge( nodeA, nodeB, lambda );
		}
	}
}

// computes the datacost of `patch` wrt to label k
int Patches::dCostB( Patch P, int k )
{
	//convert original colors to Lab, take average then reconvert avgLAB to RGB
	int maxVal = -1;
	for( int i=0 ; i<P.consistentTexColors.size() ; i++ )
	{
		cv::Vec3b f = P.consistentTexColors[i];
		int thisV = abs( (int)f[0] - k );

		if( thisV > maxVal )
			maxVal = thisV;
	}

	return maxVal;
}
// computes the datacost of `patch` wrt to label k
int Patches::dCostG( Patch P, int k )
{
	int maxVal = -1;
	for( int i=0 ; i<P.consistentTexColors.size() ; i++ )
	{
		cv::Vec3b f = P.consistentTexColors[i];
		int thisV = abs( (int)f[1] - k );

		if( thisV > maxVal )
			maxVal = thisV;
	}

	return maxVal;
}
// computes the datacost of `patch` wrt to label k
int Patches::dCostR( Patch P, int k )
{
	int maxVal = -1;
	for( int i=0 ; i<P.consistentTexColors.size() ; i++ )
	{
		cv::Vec3b f = P.consistentTexColors[i];
		int thisV = abs( (int)f[2] - k );

		if( thisV > maxVal )
			maxVal = thisV;
	}

	return maxVal;
}

// given node i and label k
int Patches::setupDataCost( Potts * p, int nLabels, int code )
{
	switch(code)
	{
	case CODE_BLUE_CHANNEL:
		for( int n=0 ; n<patchVector.size() ; n++ ) // loop over node
		{
			int * ptr = p->GetUnariesPtr(n);
			for( int k=0 ; k<nLabels ; k++ ) //loop over every possible label for this node
			{
				//ptr[k] = std::min( dCostB(patchVector[n],k), 100 ); //clipped data cost
				ptr[k] = std::min( (int)abs(patchVector[n].xavgB-k), 100 );
				//ptr[k] = std::min( abs((int)patchVector[n].getColor3dB()-k), 100 );
			}
		}
		break;


	case CODE_GREEN_CHANNEL:
		for( int n=0 ; n<patchVector.size() ; n++ ) // loop over node
		{
			int * ptr = p->GetUnariesPtr(n);
			for( int k=0 ; k<nLabels ; k++ ) //loop over every possible label for this node
			{
				//ptr[k] = std::min( dCostG(patchVector[n],k), 100 ); //clipped data cost
				ptr[k] = std::min( (int)abs(patchVector[n].xavgG-k), 100 );
				//ptr[k] = std::min( abs((int)patchVector[n].getColor3dG()-k), 100 );
			}
		}
		break;


	case CODE_RED_CHANNEL:
		for( int n=0 ; n<patchVector.size() ; n++ ) // loop over node
		{
			int * ptr = p->GetUnariesPtr(n);
			for( int k=0 ; k<nLabels ; k++ ) //loop over every possible label for this node
			{
				//ptr[k] = std::min( dCostR(patchVector[n],k), 100 ); //clipped data cost
				ptr[k] = std::min( (int)abs(patchVector[n].xavgR-k), 100 );
				//ptr[k] = std::min( abs((int)patchVector[n].getColor3dR()-k), 100 );
			}
		}
		break;

	default:
		cerr<< "[Patches] Illegal value of data cost\n\tShould be either of `CODE_RED_CHANNEL`, `CODE_GREEN_CHANNEL`, `CODE_BLUE_CHANNEL`"
				"Warning : Datacosts are not set optimization might be wrong.\n";
		break;
	}
}


void Patches::collectOptimalSolution(Potts * p, int code)
{
	switch(code)
	{
	case CODE_BLUE_CHANNEL:
		for( int i=0 ; i<patchVector.size() ; i++ )
		{
			//printf( "node[%d] = %d; avg = %d\n", i, p->GetSolution( i ), patchVector[i].avgB );
			patchVector[i].setOptimalB( p->GetSolution(i) );
		}
		break;


	case CODE_GREEN_CHANNEL:
		for( int i=0 ; i<patchVector.size() ; i++ )
		{
			//printf( "node[%d] = %d; avg = %d\n", i, p->GetSolution( i ), patchVector[i].avgG );
			patchVector[i].setOptimalG( p->GetSolution(i) );
		}
		break;


	case CODE_RED_CHANNEL:
		for( int i=0 ; i<patchVector.size() ; i++ )
		{
			//printf( "node[%d] = %d; avg = %d\n", i, p->GetSolution( i ), patchVector[i].avgR );
			patchVector[i].setOptimalR( p->GetSolution(i) );
		}
		break;
	}
}

