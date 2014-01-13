#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#ifdef HAS_GOTOOLS
	#include <GoTools/geometry/SplineSurface.h>
	#include <GoTools/trivariate/SplineVolume.h>
	#include <GoTools/geometry/ObjectHeader.h>
#endif
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/MeshRectangle.h"

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
	int it = 9;
	int dim             = 3;
	bool rat            = false;
	bool vol            = false;
	char *lrInitMesh    = NULL;
	stringstream parameters;
	parameters << " parameters: \n" \
	              "   -p1    <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	              "   -p2    <n>  polynomial order in second parametric direction\n" \
	              "   -p3    <n>  polynomial order in third parametric direction (only trivariate volumes)\n" \
	              "   -p     <n>  polynomial order in all parametric directions \n" \
	              "   -n1    <n>  number of basis functions in first parametric direction\n" \
	              "   -n2    <n>  number of basis functions in second parametric direction\n" \
	              "   -n3    <n>  number of basis functions in third parametric direction (only trivariate volumes)\n" \
	              "   -n     <n>  number of basis functions in all parametric directions\n" \
	              "   -it    <n>  number of evaluation points per element\n" \
	              "   -in:   <s>  make the LRSplineSurface <s> the initial mesh\n"\
	              "   -help       display (this) help screen\n";
	parameters << " default values\n";
	parameters << "   -p   = { " << p1 << ", " << p2 << ", " << p3 << " }\n";
	parameters << "   -n   = { " << n1 << ", " << n2 << ", " << n3 << " }\n";
	parameters << "   -it  = { " << it << " }\n";
	parameters << "   -vol = { " << ((vol)?"true":"false") << " }\n";
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0) {
			p1 = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-p2") == 0) {
			p2 = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-p3") == 0) {
			p3  = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-p") == 0) {
			p1 = p2 = p3 = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-n1") == 0) {
			n1 = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-n2") == 0) {
			n2 = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-n3") == 0) {
			n3  = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-n") == 0) {
			n1 = n2 = n3 = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-it") == 0) {
			it = atoi(argv[++i]);
		} else if(strncmp(argv[i], "-in:", 4) == 0) {
			lrInitMesh = argv[i]+4;
		} else if(strcmp(argv[i], "-vol") == 0) {
			vol = true;
		} else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << "[parameters] " << endl << parameters.str();
			exit(0);
		} else {
			cerr << "usage: " << argv[0] << "[parameters] " << endl << parameters.str();
			exit(1);
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

	LRSplineSurface *lrs;
	LRSplineVolume  *lrv;
	LRSpline        *lr;

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
			lrv = new LRSplineVolume(n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
			lr  = lrv;
		} else {
			lrs = new LRSplineSurface(n1, n2,    p1, p2,     knot_u, knot_v,         cp, dim, rat);
			lr  = lrs;
		}
	} else {
#ifdef HAS_GOTOOLS
		ifstream inputfile;
		inputfile.open(lrInitMesh);
		if(!inputfile.is_open()) {
			cerr << "Error: could not open file " << lrInitMesh << endl;
			exit(3);
		}
		Go::ObjectHeader head;
		Go::SplineSurface *ss = new Go::SplineSurface();
		Go::SplineVolume  *sv = new Go::SplineVolume();
		inputfile >> head;
		if(head.classType() == Go::Class_SplineVolume) {
			vol = true;
			inputfile >> *sv;
			lr = lrv = new LRSplineVolume(sv);
		} else if(head.classType() == Go::Class_SplineSurface) {
			vol = false;
			inputfile >> *ss;
			lr = lrs = new LRSplineSurface(ss);
		}
#endif
	}
	lr->generateIDs();


	// ---------------- Do evaluation on known elements  --------------
	vector<vector<double> > result;
	vector<double> pt;
	int maxDerivs = max(p1,p2);
	if(vol) maxDerivs = max(maxDerivs, p3);
	maxDerivs -= 2;
	char string[256];
	bool firstPrint = true;
	for(Element *el : lr->getAllElements()) {
		double u = (el->getParmin(0) + el->getParmax(0))/2.0;
		double v = (el->getParmin(1) + el->getParmax(1))/2.0;
		double w = (vol) ? (el->getParmin(2) + el->getParmax(2))/2.0 : 0.0;
		for(int i=0; i<it; i++) {
			{
			PROFILE("w/id basis");
			if(vol)
				lrv->computeBasis(u,v,w, result, 0, el->getId());
			else
				lrs->computeBasis(u,v,   result, 0, el->getId());
			}
			{
			sprintf(string, "w/id %d derivs", maxDerivs);
			PROFILE(string);
			if(vol)
				lrv->computeBasis(u,v,w, result, maxDerivs, el->getId());
			else
				lrs->computeBasis(u,v,   result, maxDerivs, el->getId());
			}
			if(firstPrint && i==0)
				cout << "w/id derivs result size: " << result.size() << " x " << result[0].size() << endl;
			{
			PROFILE("w/id point");
			if(vol)
				lrv->point(pt, u,v,w, el->getId());
			else
				lrs->point(pt, u,v,   el->getId());
			}
			{
			PROFILE("wo/id basis");
			if(vol)
				lrv->computeBasis(u,v,w, result, 0);
			else
				lrs->computeBasis(u,v,   result, 0);
			}
			{
			sprintf(string, "wo/id %d derivs", maxDerivs);
			PROFILE(string);
			if(vol)
				lrv->computeBasis(u,v,w, result, maxDerivs);
			else
				lrs->computeBasis(u,v,   result, maxDerivs);
			}
			if(firstPrint && i==0)
				cout << "wo/id derivs result size: " << result.size() << " x " << result[0].size() << endl;
			{
			PROFILE("wo/id point");
			if(vol)
				lrv->point(pt, u,v,w);
			else
				lrs->point(pt, u,v  );
			}
		}
		firstPrint = false;
	}

}
