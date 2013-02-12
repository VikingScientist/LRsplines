#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/MeshRectangle.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {

	// set default parameter values
	int p1 = 3;
	int p2 = 3;
	int p3 = 3;
	int n1 = 9;
	int n2 = 9;
	int n3 = 9;
	int dim = 4;
	bool rat = false;
	bool vol = false;
	char *inputFileName = NULL;
	string parameters(" parameters: \n" \
	                  "   -p1  <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	                  "   -p2  <n>  polynomial order in second parametric direction\n" \
	                  "   -p3  <n>  polynomial order in third parametric direction\n" \
	                  "   -n1  <n>  number of basis functions in first parametric direction\n" \
	                  "   -n2  <n>  number of basis functions in second parametric direction\n" \
	                  "   -n3  <n>  number of basis functions in third parametric direction\n" \
	                  "   -dim <n>  dimension of the controlpoints\n" \
	                  "   -vol      enforce trivariate volumetric test case\n" \
	                  "   -help     display (this) help screen\n");
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p3") == 0) {
			p3 = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n3") == 0) {
			n3 = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-dim") == 0)
			dim = atoi(argv[++i]);
		else if(strcmp(argv[i], "-vol") == 0)
			vol = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cerr << "usage: " << argv[0] << " [inputfile] [parameters]" << endl << parameters;
			exit(0);
		} else {
			if(inputFileName != NULL) {
				cerr << "usage: " << argv[0] << " [inputfile] [parameters]" << endl << parameters;
				exit(1);
			} else {
				inputFileName = argv[i];
			}
		}
	}

	// do some error testing on input
	if(n1 < p1) {
		cerr << "ERROR: n1 must be greater or equal to p1\n";
		exit(2);
	} else if(n2 < p2) {
		cerr << "ERROR: n2 must be greater or equal to p2\n";
		exit(2);
	} else if(n3 < p3) {
		cerr << "ERROR: n3 must be greater or equal to p3\n";
		exit(2);
	}

		
	// read the surface, or make one up if none specified
	LRSplineSurface *lrs = NULL;
	LRSplineVolume  *lrv = NULL;
	if(inputFileName == NULL) {
		// make a uniform integer knot vector 
		double knot_u[n1+p1];
		double knot_v[n2+p2];
		double knot_w[n3+p3];
		for(int i=0; i<p1+n1; i++)
			knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
		for(int i=0; i<p2+n2; i++)
			knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;
		for(int i=0; i<p3+n3; i++)
			knot_w[i] = (i<p3) ? 0 : (i>n3) ? n3-p3+1 : i-p3+1;
	
		// create a list of random control points (all between 0.1 and 1.1)
		int nCP  = (vol) ? n1*n2*n3 : n1*n2;
		nCP     *= (dim+rat);
		double cp[nCP];
		int k=0;
		for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
			cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0

		if(vol) {
			lrv = new LRSplineVolume(n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
			// insert a cross in the lower left corner element (should always be possible)
			lrv->insert_line(new MeshRectangle(.5,0,0,  .5,1,1) );
			lrv->insert_line(new MeshRectangle(0,.5,0,  1,.5,1) );
			lrv->insert_line(new MeshRectangle(0,0,.5,  1,1,.5) );
		} else {
			lrs = new LRSplineSurface(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
			// insert a cross in the lower left corner element (should always be possible)
			lrs->insert_const_u_edge(.5, 0, 1);
			lrs->insert_const_v_edge(.5, 0, 1);
		}


	} else {
		ifstream inputfile;
		inputfile.open(inputFileName);
		if(!inputfile.is_open()) {
			cerr << "Error: could not open file " << inputFileName << endl;
			exit(3);
		}
		lrs = new LRSplineSurface();
		inputfile >> *lrs;
	}

	// test writing to file
	ofstream lrfile;
	lrfile.open("TestReadWrite.lr");
	lrfile.precision(16);
	if(vol) lrfile << *lrv << endl;
	else    lrfile << *lrs << endl;
	
	lrfile.close();

	// test reading from file
	LRSplineSurface inputSplineSurf;
	LRSplineVolume  inputSplineVol;
	ifstream inputFile;
	inputFile.open("TestReadWrite.lr");
	if(vol) inputSplineVol.read(inputFile);
	else    inputSplineSurf.read(inputFile);
	inputFile.close();

	// write this out again and see if the are the same
	ofstream lrfile2;
	lrfile2.open("TestReadWrite2.lr");
	if(vol) lrfile2 << inputSplineVol  << endl;
	else    lrfile2 << inputSplineSurf << endl;
	lrfile2.close();

	// take a (deep) copy, screw up the original and write the copied LR spline
	// should remain unchanged if it is a proper deep copy
	ofstream lrfile3;
	ofstream lrfile4;
	lrfile3.open("TestReadWrite3.lr");
	lrfile4.open("TestReadWrite4.lr"); // this SHOULD be different from 1-3. Reg test shouldn't check against this one
	if(vol) {
		LRSplineVolume *copyVol = lrv->copy();
		lrv->getBasisfunction(0)->getknots(0)[0] = -99999;
		lrv->getElement(0)->setUmin(               -99999);
		lrv->getMeshRectangle(0)->start_[0]      = -99999;
		lrfile3 << *copyVol << endl;
		lrfile4 << *lrv << endl;
	} else {
		LRSplineSurface *copySurf = lrs->copy();
		lrs->getBasisfunction(0)->getknots(0)[0] = -99999;
		lrs->getElement(0)->setUmin(               -99999);
		(*lrs->meshlineBegin())->start_          = -99999;
		lrfile3 << *copySurf << endl;
		lrfile4 << *lrs << endl;
	}
	lrfile3.close();

}

