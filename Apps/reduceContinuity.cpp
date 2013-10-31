#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {
	if(argc < 2) {
		cout << "Usage: " << argv[0] << " <inputfile>\n";
		exit(1);
	}
	
	ifstream inputfile;
	inputfile.open(argv[1]);
	if(!inputfile.is_open()) {
		cerr << "Error: could not open file " << argv[1] << endl;
		exit(2);
	}
	LRSplineSurface lr;
	inputfile >> lr;
	int p1 = lr.order(0);
	int p2 = lr.order(1);

	for(Meshline *m : lr.getAllMeshlines()) {
		if(m->span_u_line_)
			lr.insert_const_v_edge(m->const_par_, m->start_, m->stop_, p2-1);
		else 
			lr.insert_const_u_edge(m->const_par_, m->start_, m->stop_, p1-1);
	}

	cout << lr;

}
