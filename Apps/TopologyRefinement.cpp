#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/trivariate/SplineVolume.h>
#include <GoTools/geometry/ObjectHeader.h>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/MeshRectangle.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {
#ifdef TIME_LRSPLINE
	Profiler prof(argv[0]);
#endif

	// set default parameter values
	int p1 = 3;
	int p2 = 2;
	int p3 = 4;
	int n1 = 6;
	int n2 = 5;
	int n3 = 8;
	int dim             = 3;
	bool rat            = false;
	bool dumpFile       = false;
	bool quiet          = false;
	bool vol            = false;
	char *lrInitMesh    = NULL;
	char *inputFileName = NULL;
	stringstream parameters;
	parameters << " parameters: \n" \
	              "   -p1    <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	              "   -p2    <n>  polynomial order in second parametric direction\n" \
	              "   -p3    <n>  polynomial order in third parametric direction (only trivariate volumes)\n" \
	              "   -n1    <n>  number of basis functions in first parametric direction\n" \
	              "   -n2    <n>  number of basis functions in second parametric direction\n" \
	              "   -n3    <n>  number of basis functions in third parametric direction (only trivariate volumes)\n" \
	              "   -in:   <s>  make the LRSplineSurface <s> the initial mesh\n"\
	              "   -quiet      don't dump lrspline file to stdout \n" \
	              "   -dumpfile   writes an eps- and txt-file of the LR-mesh (bivariate surfaces only)\n"\
	              "   -help       display (this) help screen\n";
	parameters << " default values\n";
	parameters << "   -p   = { " << p1 << ", " << p2 << ", " << p3 << " }\n";
	parameters << "   -n   = { " << n1 << ", " << n2 << ", " << n3 << " }\n";
	parameters << "   -vol = { " << ((vol)?"true":"false") << " }\n";
	parameters << " <refine inputfile>\n"\
	              "   inputfile describing meshline insertions.\n"\
	              "   FORMAT:\n"\
	              "     <numb. refined elements>\n"\
	              "     <insert all meshlines at once>\n"\
	              "     <x> <y>        (Surface only)\n"\
	              "     <x> <y> <z>    (Volume only)\n";
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p3") == 0) {
			p3  = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n3") == 0) {
			n3  = atoi(argv[++i]);
			vol = true;
		} else if(strncmp(argv[i], "-in:", 4) == 0)
			lrInitMesh = argv[i]+4;
		else if(strcmp(argv[i], "-vol") == 0)
			vol = true;
		else if(strcmp(argv[i], "-quiet") == 0)
			quiet = true;
		else if(strcmp(argv[i], "-dumpfile") == 0)
			dumpFile = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters.str();
			exit(0);
		} else {
			if(inputFileName != NULL) {
			cerr << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters.str();
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
	} else if(inputFileName == NULL) {
		cerr << "ERROR: Specify input file name\n";
		cerr << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters.str();
		exit(3);
	}

	SplineSurface   *ss;
	SplineVolume    *sv;
	LRSplineSurface *lrs;
	LRSplineVolume  *lrv;

	if(lrInitMesh == NULL) {
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
		int nCP = (vol) ? n1*n2*n3 : n1*n2;
		nCP    *= (dim+rat);
		double cp[nCP];
		int k=0;
		for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
			cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0
			
		// make two spline objects
		if(vol) {
			sv  = new SplineVolume  (n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
			lrv = new LRSplineVolume(n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
		} else {
			ss  = new SplineSurface  (n1, n2,    p1, p2,     knot_u, knot_v,         cp, dim, rat);
			lrs = new LRSplineSurface(n1, n2,    p1, p2,     knot_u, knot_v,         cp, dim, rat);
		}
	} else {
		ifstream inputfile;
		inputfile.open(lrInitMesh);
		if(!inputfile.is_open()) {
			cerr << "Error: could not open file " << lrInitMesh << endl;
			exit(3);
		}
		ObjectHeader head;
		ss = new SplineSurface();
		sv = new SplineVolume();
		inputfile >> head;
		if(head.classType() == Class_SplineVolume) {
			vol = true;
			inputfile >> *sv;
			lrv = new LRSplineVolume(sv);
		} else if(head.classType() == Class_SplineSurface) {
			vol = false;
			inputfile >> *ss;
			lrs = new LRSplineSurface(ss);
		}
	}

	// ---------------- Read Input file   ----------------------------
	vector<vector<double> > parPt;
	bool atOnce = false;
	ifstream inputFile;
	inputFile.open(inputFileName);
	if( inputFile.is_open() ) {
		int n;
		int parDim = (vol) ? 3 : 2;
		inputFile >> n;
		inputFile >> atOnce;
		for(int i=0; i<n; i++) {
			parPt.push_back(vector<double>(parDim));
			for(int j=0; j<parDim; j++)
				inputFile >> parPt[i][j];
		}
	}

	// ---------------- Do actual refinement   -----------------------
	if(atOnce) {
		vector<int> elIndex;
		for(uint i=0; i<parPt.size(); i++) {
			if(vol)
				elIndex.push_back(lrv->getElementContaining(parPt[i][0], parPt[i][1], parPt[i][2]) );
			else 
				elIndex.push_back(lrs->getElementContaining(parPt[i][0], parPt[i][1]) );
		}
		if(vol) lrv->refineElement(elIndex);
		else    lrs->refineElement(elIndex);
	} else {
		for(uint i=0; i<parPt.size(); i++) {
			if(vol) lrv->refineElement(lrv->getElementContaining(parPt[i][0], parPt[i][1], parPt[i][2]) );
			else    lrs->refineElement(lrs->getElementContaining(parPt[i][0], parPt[i][1]) );
		}
	}

	// ---------------- Print the full spline details  ---------------
	if(vol && !quiet)
		cout << *lrv << endl;
	else if(!vol && !quiet)
		cout << *lrs << endl;

	
	// ---------------- Print topology information   -----------------
	int nBasis;
	int nMeshlines;
	int nElements;
	if(vol) {
		nBasis     = lrv->nBasisFunctions();
		nMeshlines = lrv->nMeshRectangles();
		nElements  = lrv->nElements();
	} else {
		nBasis     = lrs->nBasisFunctions();
		nMeshlines = lrs->nMeshlines();
		nElements  = lrs->nElements();
	}
	cout << "Key LR-spline information:\n";
	cout << "  number of basis functions: " << nBasis     << endl;
	cout << "  number of mesh lines     : " << nMeshlines << endl;
	cout << "  number of elements       : " << nElements  << endl;
	
	
	// ---------------- Dump detailed debug files    -----------------
	if(dumpFile) {
		cout << endl;
		cout << "Written ";

		if(!vol) {
			ofstream meshfile;
			meshfile.open("mesh.eps");
			lrs->writePostscriptMesh(meshfile);
			meshfile.close();

			ofstream functionfile;
			functionfile.open("functions.eps");
			lrs->writePostscriptFunctionSpace(functionfile);
			functionfile.close();

			ofstream domainfile;
			domainfile.open("domain.eps");
			lrs->writePostscriptElements(domainfile, 10, 10);
			domainfile.close();

			ofstream controlmesh;
			controlmesh.open("controlmesh.eps");
			lrs->writePostscriptMeshWithControlPoints(controlmesh, 10, 10);
			controlmesh.close();
			cout << "mesh.eps, functions.eps, domain.eps, controlmesh.eps and ";
		}
		cout << "RefUnchagned.lr\n";
	
		ofstream lrfile;
		lrfile.open("RefUnchagned.lr");
		if(vol) lrfile << *lrv << endl;
		else    lrfile << *lrs << endl;
		lrfile.close();
	}

}
