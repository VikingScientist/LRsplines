#include <stdio.h>
#include <cstdlib>
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
	int n1 = 9;
	int n2 = 9;
	int dim = 4;
	bool rat = false;
	char *inputFileName = NULL;
	string parameters(" parameters: \n" \
	                  "   -p1  <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	                  "   -p2  <n>  polynomial order in second parametric direction\n" \
	                  "   -n1  <n>  number of basis functions in first parametric direction\n" \
	                  "   -n2  <n>  number of basis functions in second parametric direction\n" \
	                  "   -dim <n>  dimension of the controlpoints\n" \
	                  "   -help     display (this) help screen\n");

	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-dim") == 0)
			dim = atoi(argv[++i]);
		else if(strcmp(argv[i], "-help") == 0) {
			cerr << "usage: " << argv[0] << " [inputfile] [parameters]" << endl << parameters.c_str();
			exit(0);
		} else {
			if(inputFileName != NULL) {
				cerr << "usage: " << argv[0] << " [inputfile] [parameters]" << endl << parameters.c_str();
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
	}


	// read the surface, or make one up if none specified
	LRSplineSurface *lrs = NULL;
	if(inputFileName == NULL) {
		// make a uniform integer knot vector
		std::vector<double> knot_u(n1 + p1);
		std::vector<double> knot_v(n2 + p2);
		for(int i=0; i<p1+n1; i++)
			knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
		for(int i=0; i<p2+n2; i++)
			knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;

		// create a list of random control points (all between 0.1 and 1.1)
		int nCP  = n1*n2;
		nCP     *= (dim+rat);
		std::vector<double> cp(nCP);
		int k=0;
		for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
			cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0

		lrs = new LRSplineSurface(n1, n2, p1, p2, knot_u.begin(), knot_v.begin(), cp.begin(), dim, rat);
		// insert a cross in the lower left corner element (should always be possible)
		lrs->insert_const_u_edge(.5, 0, 1);
		lrs->insert_const_v_edge(.5, 0, 1);

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

	double knot1[] = {0,1,2,3,4};
	double knot2[] = {1,1,2,3,3};
	double cp[]    = {2,3};
	Basisfunction b1(knot1, knot2, cp, 2, 4, 4);
	vector<double> newKnot;
	b1.order_elevate(newKnot, 0);
	for(auto d : newKnot) cout << d << ", ";
	cout << endl;
	b1.order_elevate(newKnot, 1);
	for(auto d : newKnot) cout << d << ", ";
	cout << endl;



}

