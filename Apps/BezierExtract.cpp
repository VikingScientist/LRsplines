#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"

using namespace LR;
using namespace std;

int main(int argc, char **argv) {
	int el = -1;
	string filename = "";

	// read input parameters
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i],"-el") == 0) {
			el = atoi(argv[++i]);
		} else if(filename.length() > 0) {
			cout << "Usage: " << argv[0] << " [-el <n>] <inputfile>\n";
			exit(1);
		} else  {
			filename = argv[i];
		}
	}
	
	if(filename.length() == 0) {
		cout << "Usage: " << argv[0] << " [-el <n>] <inputfile>\n";
		exit(1);
	}
	
	// read lr file
	ifstream inputfile;
	inputfile.open(filename.c_str());
	if(!inputfile.is_open()) {
		cerr << "Error: could not open file " << filename << endl;
		exit(2);
	}
	LRSplineVolume lr;
	inputfile >> lr;

	if(lr.dimension() != 3) {
		cerr << "Nag on lazy Kjetil to make this test for dimensions other than 3 as well\n";
		exit(3);
	}

	int comp = 1;
	for(int i=0; i<3; i++)
		comp *= lr.order(i);
	
	vector<double> cp;
	if(el < 0) {
		for(int i=0; i<lr.nElements(); i++) {
			printf("Element #%3d\n", i);
			lr.getBezierElement(i, cp);
			int ip=0;
			for(int j=0; j<comp; j++, ip+=3)
				printf("  (%.3f, %.3f, %.3f)\n", cp[ip], cp[ip+1], cp[ip+2]);
		}
	} else {
		lr.getBezierElement(el, cp);
		printf("Element #%3d\n", el);
		int ip=0;
		for(int j=0; j<comp; j++, ip+=3)
			printf("  (%.3f, %.3f, %.3f)\n", cp[ip], cp[ip+1], cp[ip+2]);
	}

}
