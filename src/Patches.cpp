/*
 * Patches.cpp
 *
 *  Created on: 23 Dec, 2014
 *      Author: eeuser
 */

#include "Patches.h"

Patches::Patches() {
	pmvsDir = string("");

}

Patches::~Patches() {
	// TODO Auto-generated destructor stub
}

bool Patches::readPatchFile()
{
	if( pmvsDir.compare("") == 0 )
	{
		cerr<< "[Patches] PMVS directory not set\n";
		return false;
	}
	string patchFileName = pmvsDir + string("/models/option-0000.patch");

	//
	// Open File
	//
	ifstream fHandle;
	fHandle.open( patchFileName.c_str() );

	if( !fHandle.is_open() )
	{
		cerr << "[Patches] Cannot Open file :" << patchFileName << endl;
		return false;
	}

	cout<< "[Patches] Successfully opened file : "<< patchFileName << endl;


	//
	// Extract Header
	// Expecting first line to be "PATCHES"
	//
	string header;
	fHandle >> header;

	if( header.compare( "PATCHES" ) != 0 )
	{
		cerr<< "[Patches] Invalid file header\n";
		return false;
	}

	cout<< "[Patches] Header ... OK\n";


	//
	// Get total size and reserve memory
	//
	cout << "[Patches] Attempting to read rest of the file\n";
	int nSize=-1;
	fHandle >> nSize;
	if( nSize <= 0 ) { cerr<< "[Patches] Invalid nPatchs\n"; return false; }

	nPatchs = nSize;
	patchVector.reserve(nSize); //reserve memory for storing patch info

	//while (!fHandle.eof())
	for( int i=0 ; i<nPatchs ; i++)
	{
		fHandle >> header;
		if( header.compare( "PATCHS" ) == 0 )
		{
			//process this patch
			double tmp[20];
			//read 11 floats. ::: 2 vec4, consistency-score and 2 other floats
			for( int ss=0; ss<11 ; ss++ ){
				fHandle >> tmp[ss];
				//cout<< tmp[ss] << " ";
			}//cout<< "\n";
			Patch tmppatch;
			tmppatch.setPos3d(tmp[0], tmp[1], tmp[2], tmp[3]);
			tmppatch.setNormal(tmp[4], tmp[5], tmp[6], tmp[7]);
			tmppatch.setPhotometricConsistencyScore( tmp[8] );
			tmppatch.setDebug(tmp[9], tmp[10]);

			int cImIndxSize, nImIndxSize; // consistent and non-consistent image index size
			// get consistent image index
			fHandle >> cImIndxSize; int itmp;
			tmppatch.consistentTex.reserve(cImIndxSize);
			while( cImIndxSize-- > 0) {
				fHandle >> itmp;
				tmppatch.consistentTex.push_back(itmp);
			}

			// get visible but not texturally consistent im indx
			fHandle>>nImIndxSize;
			tmppatch.nonConsistentTex.reserve(nImIndxSize);
			while( nImIndxSize-->0){
				fHandle >> itmp;
				tmppatch.nonConsistentTex.push_back(itmp);
			}


			patchVector.push_back(tmppatch);
		}

	}
	cout << "[Patches] Successfully read file "<<patchFileName <<" \n";
	cout<< "[Patches] Contained "<< patchVector.size() << " patches\n";
	return true;

}

bool Patches::readImages()
{
	if( pmvsDir.compare("") == 0 )
	{
		cerr<< "[Patches] PMVS directory not set\n";
		return false;
	}
	//images are assumed to be in pmvsDir/visualize/%08d.jpg
	char imName[200];
	cv::Mat im;

	int count=0;
	for( int i=0; true ; i++)
	{
		snprintf( imName, 200, "%s/visualize/%08d.jpg", pmvsDir.c_str(),i);
		im = cv::imread( imName );
		if( im.empty() )
			break; // no more images

		imgs.push_back(im);
		count=i;
	}
	cout<< "[Patches] Read "<< count+1 << " images\n";
	return true;
}

bool Patches::showAllImages( )
{
	if( imgs.size()  <= 0 )
	{
		cerr<< "[Patches] No images to show\n";
		return false;
	}


	int cOffset=400, rOffset=400;
	char win[20];

	for( int i=0 ; i<imgs.size(); i++ )
	{
		snprintf( win, 20, "win%08d", i );
		cv::namedWindow(win,0);
		cv::resizeWindow(win, 200, 200 );

		cv::imshow( win, imgs[i] );
		cv::moveWindow( win, (i%4)*cOffset + 100, int(i/4)*rOffset );

	}
	cv::waitKey(0);

	return true;
}


bool Patches::readProjectionMatrices()
{
	if( pmvsDir.compare("") == 0 )
	{
		cerr<< "[Patches] PMVS directory not set\n";
		return false;
	}
	//projection matrices are assumed to be in pmvsDir/txt/%08d.txt
	char projFileName[200];
	cv::Mat tmpProjMat;

	int count=0;
	for( int i=0; true ; i++)
	{
		snprintf( projFileName, 200, "%s/txt/%08d.txt", pmvsDir.c_str(),i);
		ifstream projFHandle;
		projFHandle.open(projFileName);
		if( !projFHandle.is_open())
			break; //no more projection matrices to read

		// do the reading here
		string header;
		projFHandle >> header;
		if( header.compare("CONTOUR") != 0 )
		{
			cerr<< "[Patches] Invalid header in file :"<< projFileName << endl;
			cerr<< "[Patches] Not continuing any further......quiting\n";
			return false;
		}

		//if header is ok, load 12 floats
		double ftmp;
		tmpProjMat = cv::Mat::zeros(3,4,CV_64F);
		for( int f=0; f<12; f++ )
		{
			projFHandle >> ftmp;
			int rIndx = f/4; int cIndx = f%4;
			//cout<< "M["<<rIndx<<"]["<<cIndx<<"] = "<< ftmp<<endl;
			tmpProjMat.at<double>(rIndx,cIndx) = ftmp;
		}


		projFHandle.close();

		projMat.push_back(tmpProjMat.clone()); //.clone is extremely important here
		count=i;
	}
	cout<< "[Patches] Read "<< count+1 << " projection matrices\n";


	return true;


}

const string& Patches::getPMVSDir() const {
	return pmvsDir;
}

void Patches::setPMVSDir(const string& pmvsDir) {
	this->pmvsDir = pmvsDir;
}

const vector<cv::Mat>& Patches::getImgs() const {
	return imgs;
}

const Patch& Patches::getPatchAt(int ind) const {
	return patchVector[ind];
}

const vector<cv::Mat>& Patches::getProjMat() const {
	return projMat;
}

void Patches::computeAllReprojections()
{
	cout<< "[Patches] Computing Reprojections for each patch\n";
	for( int i=0 ; i<patchVector.size() ; i++ )
	{
		patchVector[i].computeReprojections(projMat);
	}
	cout<< "[Patches] .... Complete computation of reprojections\n";
}


void Patches::computeAllReprojectionsColors()
{
	cout<< "[Patches] Computing Reprojections colors for each patch\n";
	for( int i=0 ; i<patchVector.size() ; i++ )
	{
		patchVector[i].getColorsAtReprojections(imgs);
	}
	cout<< "[Patches] .... Complete computation of reprojection colors \n";
}

void Patches::showAllImagesWithReprojections()
{
	for( int i=0 ; i<patchVector.size() ; i+=100 )
	{
		if( !patchVector[i].isReProjectionsReady() )
		{
			cerr<< "[Patches] Reprojection computation not done at patch index "<< i << endl;
			continue;
		}
		cout << " +++++ patch "<<i<<" +++++\n";
		//show all images with markings
		showThisPatch( patchVector[i] ); //displays all images with markings
		cv::waitKey(0);
	}

}

void Patches::write4Colorization()
{


	// write 3-channel gray with scribbles.
	//first check for all patch that they have a valid color, cord
	for( int i=0 ; i<patchVector.size() ; i++ )
	{
		if( !patchVector[i].isReProjectionsReady() || !patchVector[i].isReProjectionsColorsReady() )
		{
			cerr<< "[Patches] For the patch#"<< i << " reprojections and colors not computed...Quiting this function\n";
			return;
		}
	}
	cout<< "[Patches] All patches have a valid color\n";


	// create 3-channel grays for scibbles
	vector<cv::Mat> gray3Channel;
	vector<cv::Mat> chn;
	chn.reserve(3);
	for( int i=0 ; i<imgs.size() ; i++ )
	{
		cv::Mat tmpMat, tmpMat2;
		cv::cvtColor( imgs[i], tmpMat, CV_BGR2GRAY );
		//tmpMat = cv::Mat::zeros(imgs[i].rows, imgs[i].rows, CV_8UC1 );

		chn.clear();
		chn.push_back(tmpMat);
		chn.push_back(tmpMat);
		chn.push_back(tmpMat);

		cv::merge(chn, tmpMat2);

		gray3Channel.push_back( tmpMat2.clone() );
	}

	//write all images in gray
	char outGrayFileName[200];
	cout<< "[Patches] Writing all gray images";
	for( int i=0 ; i<imgs.size() ; i++ )
	{
		snprintf( outGrayFileName, 200, "%s/visualize/gray_%08d.png", pmvsDir.c_str(), i );

		cv::imwrite( outGrayFileName, gray3Channel[i] );
	}
	cout<< " .... done!\n";


	// insert scribbles from each patch (every 100th patch actually)
	cout<< "[Patches] Insert scribbles";
	for( int i=0 ; i<patchVector.size() ; i++ )
	{
		patchVector[i].insertScribbles( gray3Channel );
	}
	cout<< " .... done!\n";


	//write scibbles
	cout<< "[Patches] Write all scribble images";
	for( int i=0 ; i<gray3Channel.size() ; i++ )
	{
		snprintf( outGrayFileName, 200, "%s/visualize/scribbles_%08d.png", pmvsDir.c_str(), i );

		cv::imwrite( outGrayFileName, gray3Channel[i] );
	}



	// write pset (
//	ply
//	format ascii 1.0
//	element vertex 241099
//	property float x
//	property float y
//	property float z
//	property float nx
//	property float ny
//	property float nz
//	property uchar diffuse_red
//	property uchar diffuse_green
//	property uchar diffuse_blue
//	end_header


	char outply[200];
	snprintf( outply, 200, "%s/models/xx.ply", pmvsDir.c_str() );
	cout << "Now writing : "<< outply << endl;
	ofstream ff(outply,std::ios::trunc);
	if( !(ff.is_open()) )
	{
		cerr<< "[Patches/write4Colorization] Cannot open file. "<< outply << endl;
		return;
	}

	ff << "ply\n";
	ff<<"format ascii 1.0\n";
	ff<<"element vertex "<<patchVector.size()<<endl;
	ff<<"property float x\n";
	ff<<"property float y\n";
	ff<<"property float z\n";
	ff<<"property float nx\n";
	ff<<"property float ny\n";
	ff<<"property float nz\n";
	ff<<"property uchar diffuse_red\n";
	ff<<"property uchar diffuse_green\n";
	ff<<"property uchar diffuse_blue\n";
	ff<<"end_header\n";

	for( int i=0 ; i<patchVector.size() ; i++ )
	{
		cv::Mat pos = patchVector[i].getPos3d();
		ff<<pos.at<double>(0)<<" "<<pos.at<double>(1)<<" "<<pos.at<double>(2)<<" ";
		cv::Mat nor = patchVector[i].getNormal();
		ff<<nor.at<double>(0)<<" "<<nor.at<double>(1)<<" "<<nor.at<double>(2)<<" ";
		ff<<patchVector[i].getOptimalR() << " "<< patchVector[i].getOptimalG()
						<< " " << patchVector[i].getOptimalB() << endl;
		//ff<<(int)patchVector[i].xavgR << " "<< (int)patchVector[i].xavgG
		//						<< " " << (int)patchVector[i].xavgB  << endl;

		//ff<< patchVector[i].getColor3dR() << " " << patchVector[i].getColor3dG()
		//		<< " " << patchVector[i].getColor3dB() << endl;
	}
	ff.close();



	cout<< " .... done!\n";


}

// read file (ply) option-colored.ply
void Patches::readPatchColors() {
	string fName = pmvsDir + string("/models/option-0000.ply");
	fstream fp(fName.c_str(), ios::in);

	if( !(fp.is_open()) )
	{
		cerr<< "[Patches::readPatchColors] Cannot open 3d color file : "<< fName;
		return;
	}

	cout<< "[Patches] Reading colored point cloud\n";
	string header;
	for( int l=0; l<6 ; l++ )
		fp>>header; //ply .... element vertex


	int totalLines=0;
	fp>>totalLines;  // 250076
	cout<< "[Patches] " << totalLines << " 3d points to read"<< endl;

	for( int l=0; l<28 ; l++ )
			fp>>header; //ply .... element vertex

	for( int i=0 ; i<totalLines ; i++ )
	{
		float tmp;
		//ignore 1st 6 numbers
		for( int x=0 ; x<6 ; x++ )
			fp>>tmp;

		int rx, gx, bx;
		fp>>rx;
		fp>>gx;
		fp>>bx;

		patchVector[i].setColor3dR(rx);
		patchVector[i].setColor3dG(gx);
		patchVector[i].setColor3dB(bx);
	}
}

void Patches::showThisPatch(Patch ch)
{
	if( imgs.size()  <= 0 )
	{
		cerr<< "[Patches] No images to show\n";
		return ;
	}

	if( ch.isReProjectionsReady() == false )
	{
		cerr<<"Reprojections are not read. Call computeAllReprojections() before calling this function\n"
				"\t not showing.....";
		return;
	}


	int cOffset=400, rOffset=400;
	char win[20];

	//assuming the win are already created
	for( int i=0 ; i<imgs.size(); i++ )
	{
		snprintf( win, 20, "win%08d", i );

		cv::imshow( win, imgs[i] );

	}


	//reshow images which has reprojections
	for( int i=0 ; i<ch.consistentTex.size() ; i++ )
	{
		int imIndx = ch.consistentTex[i];
		snprintf( win, 20, "win%08d", imIndx );
		cout<< "-----------\n";
		cv::Mat cloned = imgs[imIndx].clone();
		cv::Vec3b thisColor = ch.consistentTexColors[i];
		cv::circle(cloned, ch.consistentTexCords[i], 10, cv::Scalar((double)thisColor[0],(double)thisColor[1],(double)thisColor[2]), -1 );
		//cout<< "consistent cord : "<<imIndx<< " "<< hex<<ch.consistentTexCords[i] << endl;
		cout<< "consistent color : "<<imIndx<< " "<< hex<<ch.consistentTexColors[i] << dec<< ch.consistentTexColors[i] << endl;

		if( ch.nonConsistentTex.size() > i ){ // also draw non consistent ones, assuming non-consistent ones are smaller in number
			cv::circle(cloned, ch.nonConsistentTexCords[i], 8, cv::Scalar(0,0,255), -1 );
			//cout<< "non consistent cord : "<< imIndx<< ch.nonConsistentTexCords[i] << endl;
			cout<< "non consistent color : "<<imIndx<< hex<<ch.nonConsistentTexColors[i] << dec<<ch.nonConsistentTexColors[i]<< endl;
		}


		cv::imshow( win, cloned );
	}
}


// using http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void Patches::rgb2lab( float R, float G, float B, float & l_s, float &a_s, float &b_s )
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
void Patches::lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B )
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

