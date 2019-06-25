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
	double colors[][3] = {  {     0, 0.4470, 0.7410},
				{0.8500, 0.3250, 0.0980},
				{0.9290, 0.6940, 0.1250},
				{0.4940, 0.1840, 0.5560},
				{0.4660, 0.6740, 0.1880},
				{0.3010, 0.7450, 0.9330},
				{0.6350, 0.0780, 0.1840},
				{234.0/255, 161.0/255, 230.0/255},
				{201.0/255, 216.0/255, 32.0/255},
				{127.0/255, 104.0/255, 35.0/255},
				{160.0/255, 160.0/255, 160.0/255},
                };
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

	int pmax = lr.max_order(0);
	vector<vector<int> > myFunc(pmax+1);
	for(auto el : lr.getAllElements()) {
		myFunc[el->order(0)].push_back(el->getId());
	}

	ofstream meshfile;
	meshfile.open("mesh.eps");
	for(int p=0; p<=pmax; p++) {
		if(myFunc[p].size() == 0) continue;
		lr.setElementColor(colors[p][0], colors[p][1], colors[p][2]);
		lr.writePostscriptMesh(meshfile, false, &myFunc[p]);
	}
	lr.writePostscriptMesh(meshfile);
	meshfile.close();

	cout << "Written mesh to mesh.eps\n";
	cout << "  p=0,  blue\n";
	cout << "  p=1,  red\n";
	cout << "  p=2,  orange\n";
	cout << "  p=3,  purple\n";
	cout << "  p=4,  green\n";
	cout << "  p=5,  cyan\n";
	cout << "  p=6,  maroon\n";
	cout << "  p=7,  pink\n";
	cout << "  p=8,  lime\n";
	cout << "  p=9,  brown\n";
	cout << "  p=10, gray\n";
}
