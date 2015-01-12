/*
 * main.cpp
 *
 *  Created on: 23 Dec, 2014
 *      Author: eeuser
 */

#include <iostream>
#include <opencv2/opencv.hpp>
using namespace std;

#include "Patches.h"
#include "LBPDenoise.h"


#include "Potts.h"

void testPotts();


int main( int argc, char ** argv )
{
	/*
		LoopyBPDenoising lbp( "noisy.png" );
		lbp.denoise();
		//cv::imshow( "win", lbp.getIm() );
		//cv::waitKey(0);
		*/


	//string pmvsDir = string( "/home/eeuser/Reconstruction_Kit/a_imData/darkcorner/pmvs/");
	//string pmvsDir = string( "/home/eeuser/Reconstruction_Kit/a_imData/et/");
	string pmvsDir = string( "/home/eeuser/Reconstruction_Kit/a_imData/floor5_full/pmvs/");
	//string pmvsDir = string( "/home/eeuser/Reconstruction_Kit/a_imData/coke/pmvs/");
	Patches p;
	p.setPMVSDir(pmvsDir);
	p.readPatchFile( );
	p.readImages( );
	p.readProjectionMatrices( );



	//p.readPatchColors();
	p.computeAllReprojections();
	p.computeAllReprojectionsColors();


	//p.showAllImages();
	//p.showAllImagesWithReprojections();

	// ANN
	//testANN(p);
	p.computeSpacialNearestNeighbourGraph();


	// Graph Optimization
	//testPotts();
	p.synchronizeColors();


	p.write4Colorization();


	/*
	//testing average color
	cv::Mat im, im_out;
	im = cv::imread( "Lenna.png");
	cv::cvtColor(im, im_out, CV_BGR2Lab);
	cout<< im_out.type() << endl;
	cout<< im_out.channels() << endl;

	cv::imshow("x", im);

	cv::waitKey(0);
*/



	/*
	//testing RGB2Lab & Lab2RGB
	float R,G,B, L, a, b;
	//R=0; G=0 ; B=180;
	L=26; a=90; b= 10.;

	Patches::lab2rgb(L,a,b,R,G,B);
	cout<< round(R) << ", "<< round(G) << ", "<< round(B) << endl;

	Patches::rgb2lab(R,G,B, L, a, b );
	cout<< round(L) << ", "<< round(a) << ", "<< round(b) << endl;
*/





	cout<< "All Done....!\n";
	return 0;
}

// given node i and label k
int func( int i, int k )
{
    if( i==0 )
        return abs(200-k);
    else if( i==1 )
        return abs( 220-k );
    else if( i== 2 )
        return abs( 205-k );
    else if( i== 3 )
            return abs( 201-k );
}

void testPotts()
{
	cout<< "\nTesting potts implementation\n";

	Potts * p = new Potts( 4, 5, 256, NULL );
    p->AddEdge( 0, 1, 1 );
    p->AddEdge( 1, 2, 5 );
    p->AddEdge( 2, 3, 5 );
    p->AddEdge( 0, 3, 1 );

    int nodeCost[4] = { 200, 220, 205, 201 };
    for( int n=0 ; n<4 ; n++ ) // loop over node
    {
		int * ptr = p->GetUnariesPtr(n);
    	for( int k=0 ; k<256 ; k++ ) //loop over every possible label for this node
    	{
    		ptr[k] = abs( nodeCost[n] - k );
    	}
    }

    int c_star = p->Solve( 1, false );

    printf( "Total Enegergy = %d\n", c_star );
    for( int i=0 ; i<4 ; i++ )
        printf( "node[%d] = %d\n", i, p->GetSolution( i ) );

    return ;

}
