#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
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
			cerr << "usage: " << argv[0] << " [parameters]" << endl << parameters.c_str();
			exit(0);
		} else {
			cerr << "usage: " << argv[0] << " [parameters]" << endl << parameters.c_str();
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


	// read the surface, or make one up if none specified
	LRSpline        *lr  = NULL;
	LRSplineSurface *lrs = NULL;
	LRSplineVolume  *lrv = NULL;
	// make a uniform integer knot vector
	std::vector<double> knot_u(n1 + p1);
	std::vector<double> knot_v(n2 + p2);
	std::vector<double> knot_w(n3 + p3);
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;
	for(int i=0; i<p3+n3; i++)
		knot_w[i] = (i<p3) ? 0 : (i>n3) ? n3-p3+1 : i-p3+1;

	// create a list of random control points (all between 0.1 and 1.1)
	int nCP  = (vol) ? n1*n2*n3 : n1*n2;
	nCP     *= (dim+rat);
	std::vector<double> cp(nCP);
	int k=0;
	for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
		cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0

	if(vol) {
		lr = lrv = new LRSplineVolume(n1, n2, n3, p1, p2, p3, knot_u.begin(), knot_v.begin(), knot_w.begin(), cp.begin(), dim, rat);
	} else {
		lr = lrs = new LRSplineSurface(n1, n2, p1, p2, knot_u.begin(), knot_v.begin(), cp.begin(), dim, rat);
	}
	// refine the lower-left corner 4 times
	for(int i=0; i<4; i++) {
		lr->generateIDs();
		vector<double> corner(lr->nVariate(), 0.00001);
		Element* el = lr->getElement(lr->getElementContaining(corner));
		vector<int> refI(0);
		for(auto b : el->support()) refI.push_back(b->getId());
		lr->refineBasisFunction(refI);
	}
	lr->generateIDs();

	bool allOK = true;
	for(auto el : lr->getAllElements()) {
		double du = el->umax() - el->umin();
		int expectedLevel = round(log2(1/du));
		cout << "Element #" << el->getId() << ": Length=" << du << " level " << el->getLevel(0);
		for(int i=1; i<lr->nVariate(); i++) cout << ", " << el->getLevel(i);
		cout << " (expects lvl " << expectedLevel << ")";
		bool ok = true;
		for(int i=0; i<lr->nVariate(); i++)
			ok &= expectedLevel == el->getLevel(i);
		if(ok) cout << " - OK" << endl;
		else   cout << " - ASSERTION FAIL" << endl;
		allOK &= ok;
	}
	cout << "====================================================" << endl;
	if(allOK) cout << "All assertions passed" << endl;
	else      cout << "FAILED TEST"           << endl;
	cout << "----------------------------------------------------" << endl;
}

