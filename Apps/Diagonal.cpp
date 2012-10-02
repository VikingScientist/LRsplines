
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <GoTools/geometry/SplineSurface.h>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {

	// set default parameter values
	int p         = 3;
	int N         = 6;
	int m         = 1;
	int scheme    = 0;
	bool dumpfile = false;
	string parameters(" parameters: \n" \
	                  "   -p      <n> polynomial DEGREE (order-1) of the basis\n" \
	                  "   -n      <n> number of iterations\n" \
	                  "   -m      <n> knot multiplicity\n" \
	                  "   -dumpfile   create .eps and .m files of the meshes\n" \
	                  "   -scheme <n> refinement scheme (0=FULLSPAN, 1=MINSPAN, 2=STRUCT)\n" \
	                  "   -help    display (this) help information\n");
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p") == 0)
			p = atoi(argv[++i]);
		else if(strcmp(argv[i], "-scheme") == 0)
			scheme = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n") == 0)
			N = atoi(argv[++i]);
		else if(strcmp(argv[i], "-m") == 0)
			m = atoi(argv[++i]);
		else if(strcmp(argv[i], "-dumpfile") == 0)
			dumpfile = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << endl << parameters;
			exit(0);
		} else {
			cout << "usage: " << argv[0] << endl << parameters;
			exit(1);
		}
	}

	// do some error testing on input
	if(m > p )
		m = p;
	if(scheme > 2)
		scheme = 2;

	// clear previous result from file
	if(dumpfile) {
		ofstream clear_out;
		clear_out.open("dof.m");
		clear_out.close();
		clear_out.open("elements.m");
		clear_out.close();
	}



	// make a uniform integer knot vector
	double knot_u[] = {0,0,1,1};
	double knot_v[] = {0,0,1,1};
	double cp[]     = {0,0,
	                   1,0,
	                   0,1,
	                   1,1};
		
	// make two identical surfaces
	SplineSurface   ss(2, 2, 2, 2, knot_u, knot_v, cp, 2, false);
	ss.raiseOrder(p-1, p-1);
	LRSplineSurface lr(&ss);

	// setup refinement parameters
	lr.setRefMultiplicity(m);
	if(scheme == 0)
		lr.setRefStrat(LR_SAFE);
	else if(scheme == 1)
		lr.setRefStrat(LR_MINSPAN);
	else if(scheme == 2)
		lr.setRefStrat(LR_ISOTROPIC_FUNC);
	
	// for all iterations
	for(int n=0; n<N; n++) {
		lr.generateIDs();
		vector<int> indices;

		if(scheme < 2) {
			lr.getDiagonalElements(indices);
			lr.refineElement(indices);
		} else if(scheme == 2) {
			lr.getDiagonalBasisfunctions(indices);
			lr.refineBasisFunction(indices);
		}

		// draw result files (with next refinement-step-diagonal shaded)
		if(dumpfile) {
			vector<int> diagElms, diagFuncs;
			lr.getDiagonalElements(diagElms);
			lr.getDiagonalBasisfunctions(diagFuncs);
			char filename[128];

			sprintf(filename, "index_%02d.eps", n+1);
			ofstream out;
			out.open(filename);
			if(scheme < 2)
				lr.writePostscriptMesh(out, true, &diagElms);
			else
				lr.writePostscriptMesh(out);
			out.close();

			sprintf(filename, "parameter_%02d.eps", n+1);
			out.open(filename);
			if(scheme < 2)
				lr.writePostscriptElements(out, 2,2, true, &diagElms);
			else 
				lr.writePostscriptElements(out, 2,2, true);
			out.close();

			sprintf(filename, "func_%02d.eps", n+1);
			out.open(filename);
			if(scheme == 2)
				lr.writePostscriptFunctionSpace(out, &diagFuncs);
			else
				lr.writePostscriptFunctionSpace(out);
			out.close();
			
			out.open("dof.m", ios::app);
			out << lr.nBasisFunctions() << endl;
			out.close();

			out.open("elements.m", ios::app);
			out << lr.nElements() << endl;
			out.close();
		}
	}
	
	vector<int> overloadedBasis;
	vector<int> overloadedElements;
	vector<int> multipleOverloadedElements;

	// harvest some statistics and display these results
	lr.generateIDs();
	vector<Basisfunction*>::iterator bit;
	double avgBasisToElement = 0;
	double avgBasisToLine    = 0;
	int maxBasisToElement    = -1;
	int minBasisToElement    = 9999999;
	int maxBasisToLine       = -1;
	int minBasisToLine       = 9999999;
	int nOverloadedElms      = 0;
	int nOverloadedBasis     = 0;
	for(bit=lr.basisBegin(); bit!=lr.basisEnd(); bit++) {
		int nE = (*bit)->nSupportedElements();
		int nL = (*bit)->nPartialLines();
		maxBasisToElement = (maxBasisToElement > nE) ? maxBasisToElement : nE;
		minBasisToElement = (minBasisToElement < nE) ? minBasisToElement : nE;
		avgBasisToElement += nE;
		maxBasisToLine = (maxBasisToLine > nL) ? maxBasisToLine : nL;
		minBasisToLine = (minBasisToLine < nL) ? minBasisToLine : nL;
		avgBasisToLine += nL;
		if((*bit)->isOverloaded()) {
			nOverloadedBasis++;
			overloadedBasis.push_back((*bit)->getId());
		}
	}
	avgBasisToElement /= lr.nBasisFunctions();
	avgBasisToLine    /= lr.nBasisFunctions();

	vector<Element*>::iterator eit;
	double avgElementToBasis       = 0;
	double avgSquareElementToBasis = 0;
	int maxElementToBasis          = -1;
	int minElementToBasis          = 9999999;
	for(eit=lr.elementBegin(); eit!=lr.elementEnd(); eit++) {
		int nB = (*eit)->nBasisFunctions();
		maxElementToBasis       = (maxElementToBasis > nB) ? maxElementToBasis : nB;
		minElementToBasis       = (minElementToBasis < nB) ? minElementToBasis : nB;
		avgElementToBasis       += nB;
		avgSquareElementToBasis += nB*nB;
		if((*eit)->isOverloaded()) {
			nOverloadedElms++;
			overloadedElements.push_back((*eit)->getId());
			// int nCount = (*eit)->overloadedBasisCount();
			// if(nCount >= 2)
				// multipleOverloadedElements.push_back((*eit)->getId());
		}
	}
	avgElementToBasis /= lr.nElements();
	avgSquareElementToBasis /= lr.nElements();
	
	cout << "Some statistics: " << endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Number of basisfunctions: " << lr.nBasisFunctions()  << endl;
	cout << "Number of elements      : " << lr.nElements()        << endl;
	cout << "Number of meshlines     : " << lr.nMeshlines()       << endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Min number of Basisfuntion -> Element: " << minBasisToElement << endl;
	cout << "Max number of Basisfuntion -> Element: " << maxBasisToElement << endl;
	cout << "Avg number of Basisfuntion -> Element: " << avgBasisToElement << endl;
	cout << endl;
	cout << "Min number of Basisfuntion -> Meshline: " << minBasisToLine << endl;
	cout << "Max number of Basisfuntion -> Meshline: " << maxBasisToLine << endl;
	cout << "Avg number of Basisfuntion -> Meshline: " << avgBasisToLine << endl;
	cout << endl;
	cout << "Min number of        Element -> Basisfunction: " << minElementToBasis  << endl;
	cout << "Max number of        Element -> Basisfunction: " << maxElementToBasis  << endl;
	cout << "Avg number of        Element -> Basisfunction: " << avgElementToBasis  << endl;
	cout << "Avg square number of Element -> Basisfunction: " << avgSquareElementToBasis
	     << " (" << sqrt(avgSquareElementToBasis) << ")" << endl;
	cout << endl;
	cout << "Number of overloaded Basisfunctions : " << nOverloadedBasis  << endl;
	cout << "Number of overloaded Elements       : " << nOverloadedElms   << endl;
	cout << "-------------------------------------------------------------" << endl;
	if(lr.nBasisFunctions() < 1300) {
	cout << "Is linearly independent : " << ((lr.isLinearIndepByMappingMatrix(false) )? "True":"False") << endl;
	cout << "-------------------------------------------------------------" << endl;
	}


	if(dumpfile) {
		ofstream out;
	
		out.open("overloaded_elements.eps");
		lr.writePostscriptMesh(out, false, &overloadedElements);
		lr.setElementColor(0.9, 0.3, 0.15);
		lr.writePostscriptMesh(out, true,  &multipleOverloadedElements);
		out.close();
	
		out.open("overloaded_functions.eps");
		lr.writePostscriptFunctionSpace(out, &overloadedBasis);
		out.close();
	}
}

