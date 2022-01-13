#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#ifdef HAS_GOTOOLS
	#include <GoTools/geometry/SplineSurface.h>
	#include <GoTools/geometry/ObjectHeader.h>
#endif

#include "LRSpline/LRSpline.h"
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"

using namespace LR;
using namespace std;

int main(int argc, char **argv) {
	vector<int> el(1,0);
	int m      = 1;
	int scheme = 2;
	bool dumpfile = false;
	string filename("");
	string outfile("result.lr");
	stringstream parameters;
	parameters << " parameters: \n"
	           << "   -el     <n[]> comma-seperated list of element (or basisfunction) indices \n"
	           << "   -m       <n>  knot multiplicity\n"
	           << "   -out     <s>  output file name\n"
	           << "   -scheme  <n>  refinement scheme (0=FULLSPAN, 1=MINSPAN, 2=STRUCT)\n"
	           << "   -dumpfile     create .eps and .m files of the meshes (Surfaces only)\n"
	           << "   -help         displays (this) help screen\n"
	           << " default values:\n"
	           << "   el     = {" << el[0]   << "}\n"
	           << "   m      = "  << m       <<  "\n"
	           << "   scheme = "  << scheme  <<  "\n"
	           << "   out    = "  << outfile <<  "\n";

	// read input arguments
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-scheme") == 0)
			scheme = atoi(argv[++i]);
		else if(strcmp(argv[i], "-m") == 0)
			m = atoi(argv[++i]);
		else if(strcmp(argv[i], "-out") == 0)
			outfile = argv[++i];
		else if(strcmp(argv[i], "-el") == 0) {
			el.clear();
			char *tok = strtok(argv[++i], ",");
			while(tok != NULL) {
				el.push_back(atoi(tok));
				tok = strtok(NULL, ",");
			}
		} else if(strcmp(argv[i], "-dumpfile") == 0)
			dumpfile = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << " [parameters] " << endl << parameters.str() << endl;
			exit(0);
		} else if(filename.length() > 0) {
			cout << "usage: " << argv[0] << " [parameters] " << endl << parameters.str() << endl;
			exit(1);
		} else {
			filename = argv[i];
		}
	}

	// error test input
	if(filename.length() == 0) {
		cout << "usage: " << argv[0] << " [parameters] " << endl << parameters.str() << endl;
		exit(1);
	}
	if(scheme < 0)
		scheme = 0;
	if(scheme > 2)
		scheme = 2;

	// read input file
	ifstream inputfile;
	inputfile.open(filename.c_str());
	if(!inputfile.is_open()) {
		cerr << "Error: could not open file " << filename << endl;
		exit(2);
	}
	LRSpline        *lr;
	LRSplineSurface *lrs = NULL;
	LRSplineVolume  *lrv = NULL;
  int p = 0;
	char buffer[512];
	inputfile.getline(buffer, 512); // peek the first line to figure out if it's an LRSpline or a GoTools spline
	inputfile.seekg(ios_base::beg);
	if(strncmp(buffer, "# LRSPLINE VOLUME",17)==0) {
		lr = lrv = new LRSplineVolume();
		inputfile >> *lrv;
		p = std::min(std::min(lrv->min_order(0), lrv->min_order(1)), lrv->min_order(2));
	} else if(strncmp(buffer, "# LRSPLINE",10)==0) {
		lr = lrs = new LRSplineSurface();
		inputfile >> *lrs;
		p = std::min(lrs->min_order(0), lrs->min_order(1));
	} else {
#ifdef HAS_GOTOOLS
		Go::ObjectHeader   head;
		inputfile >> head;
		if(head.classType() == Go::Class_SplineVolume) {
			Go::SplineVolume     sv;
			inputfile >> sv;
			lr = lrv = new LRSplineVolume(&sv);
		} else if(head.classType() == Go::Class_SplineSurface) {
			Go::SplineSurface    ss;
			inputfile >> ss;
			lr = lrs = new LRSplineSurface(&ss);
		}  else {
			std::cerr << "Unsupported GoTools object\n";
			exit(3);
		}
#endif
	}

	// setup refinement parameters
	lr->setRefContinuity(p-m-1);
	if(scheme == 0)
		lr->setRefStrat(LR_FULLSPAN);
	else if(scheme == 1)
		lr->setRefStrat(LR_MINSPAN);
	else if(scheme == 2)
		lr->setRefStrat(LR_STRUCTURED_MESH);

	// do actual refinement
	if(scheme == 2)
		lr->refineBasisFunction(el);
	else
		lr->refineElement(el);

	// write result to file
	ofstream refinedFile;
	refinedFile.open(outfile.c_str());
	refinedFile << *lr;
	refinedFile.close();
	cout << "Written refined LR-spline to " << outfile << endl;

	// dump debug eps files
	if(dumpfile && lrs != NULL) {
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

		cout << endl;
		cout << "Written mesh to mesh.eps, functions.eps, domain.eps and controlmesh.eps \n";
	}
}
