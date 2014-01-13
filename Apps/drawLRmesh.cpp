#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"

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

	ofstream meshfile;
	meshfile.open("mesh.eps");
	lr.writePostscriptMesh(meshfile);
	meshfile.close();

	ofstream functionfile;
	functionfile.open("functions.eps");
	lr.writePostscriptFunctionSpace(functionfile);
	functionfile.close();

	ofstream domainfile;
	domainfile.open("domain.eps");
	lr.writePostscriptElements(domainfile, 10, 10);
	domainfile.close();

	ofstream controlmesh;
	controlmesh.open("controlmesh.eps");
	lr.writePostscriptMeshWithControlPoints(controlmesh, 10, 10);
	controlmesh.close();

	cout << endl;
	cout << "Written mesh to mesh.eps, functions.eps, domain.eps and controlmesh.eps \n";
}
