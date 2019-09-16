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
	double colors[][3] = {
	            {0.8600, 0.8600, 0.8600},
	            {     0, 0.4470, 0.7410}, // default matlab colors from here (p=1)
	            {0.8500, 0.3250, 0.0980},
	            {0.9290, 0.6940, 0.1250},
	            {0.4940, 0.1840, 0.5560},
	            {0.4660, 0.6740, 0.1880},
	            {0.3010, 0.7450, 0.9330},
	            {0.6350, 0.0780, 0.1840}, // and down to here (p=7)
	            {234.0/255, 161.0/255, 230.0/255},
	            {201.0/255, 216.0/255, 32.0/255},
	            {127.0/255, 104.0/255, 35.0/255},
	            {140.0/255, 140.0/255, 140.0/255},
	             };
	if(argc < 2) {
		cout << "Usage: " << argv[0] << " <inputfile>\n";
		exit(1);
	}
    bool silent = false;

	ifstream inputfile;
	inputfile.open(argv[1]);
	if(!inputfile.is_open()) {
		cerr << "Error: could not open file " << argv[1] << endl;
		exit(2);
	}
	if(argc > 2 && strcmp(argv[2], "--silent")==0)
		silent = true;

	LRSplineSurface lr;
	inputfile >> lr;

	int pmax = lr.max_order(0);
	vector<vector<int> > myFunc(pmax+1);
	for(auto el : lr.getAllElements()) {
		myFunc[el->order(0)].push_back(el->getId());
	}

	ofstream meshfile;
	meshfile.open("mesh.eps");
	ofstream elementfile;
	elementfile.open("elements.eps");
	for(int p=1; p<=pmax; p++) {
		if(myFunc[p].size() == 0) continue;
		lr.setElementColor(colors[p-1][0], colors[p-1][1], colors[p-1][2]);
		lr.writePostscriptMesh(    meshfile,        false, &myFunc[p]);
		lr.writePostscriptElements(elementfile, 3,3,false, &myFunc[p]);
	}
	lr.writePostscriptMesh(meshfile);
	lr.writePostscriptElements(elementfile,3,3,true);
	meshfile.close();
	elementfile.close();

    if(!silent) {
		cout << "Written computational mesh to mesh.eps\n";
		cout << "Written domain to elements.eps\n";
		cout << "  p=0,  light gray\n";
		cout << "  p=1,  blue\n";
		cout << "  p=2,  red\n";
		cout << "  p=3,  orange\n";
		cout << "  p=4,  purple\n";
		cout << "  p=5,  green\n";
		cout << "  p=6,  cyan\n";
		cout << "  p=7,  maroon\n";
		cout << "  p=8,  pink\n";
		cout << "  p=9,  yellow\n";
		cout << "  p=10, brown\n";
		cout << "  p=11, dark gray\n";
	}
}
