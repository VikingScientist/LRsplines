#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {
	Profiler prof(argv[0]);

	// set default parameter values
	int p1 = 3;
	int p2 = 3;
	int n1 = 9;
	int n2 = 9;
	int dim = 4;
	bool rat = false;
	string parameters(" parameters: \n" \
	                  "   -p1  <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	                  "   -p2  <n>  polynomial order in second parametric direction\n" \
	                  "   -n1  <n>  number of basis functions in first parametric direction\n" \
	                  "   -n2  <n>  number of basis functions in second parametric direction\n" \
	                  "   -dim <n>  dimension of the controlpoints\n" \
	                  "   -help     display (this) help screen");
	
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
			cerr << "usage: " << argv[0] << endl << parameters;
			exit(0);
		} else {
			cerr << "usage: " << argv[0] << endl << parameters;
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
	}

	// make a uniform integer knot vector
	double knot_u[n1+p1];
	double knot_v[n2+p2];
	double cp[(dim+rat)*n1*n2];
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;

	// create a list of random control points (all between 0 and 12)
	int k=0;
	for(int j=0; j<n2; j++) 
		for(int i=0; i<n1; i++) 
			for(int d=0; d<dim+rat; d++)
				cp[k++] = ((i*2+j*3+d*5) % 13); 
		
	// make two identical surfaces
	LRSplineSurface lr(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);

	// insert a cross in the lower left corner element (should always be possible)
	lr.insert_const_u_edge(.5, 0, 1);
	lr.insert_const_v_edge(.5, 0, 1);

	// test writing to file
	ofstream lrfile;
	lrfile.open("TestReadWrite.lr");
	lrfile << lr << endl;
	lrfile.close();

	// test reading from file
	LRSplineSurface inputSpline;
	ifstream inputFile;
	inputFile.open("TestReadWrite.lr");
	inputSpline.read(inputFile);
	inputFile.close();

	// write this out again and see if the are the same
	ofstream lrfile2;
	lrfile2.open("TestReadWrite2.lr");
	lrfile2 << inputSpline << endl;
	lrfile2.close();
}

