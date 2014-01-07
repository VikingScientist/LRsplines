
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
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
	bool vol       = false;
	string parameters(" parameters: \n" \
	                  "   -p      <n> polynomial DEGREE (order-1) of the basis\n" \
	                  "   -n      <n> number of iterations\n" \
	                  "   -m      <n> knot multiplicity\n" \
	                  "   -dumpfile   create .eps and .m files of the meshes\n" \
	                  "   -vol        enforce a volumetric test case\n"\
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
		else if(strcmp(argv[i], "-vol") == 0)
			vol = true;
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



	// make two identical splines
	LRSplineSurface *lrs;
	LRSplineVolume  *lrv;
	if(vol) {
		lrv = new LRSplineVolume(p+1, p+1, p+1, p+1, p+1, p+1);
	} else {
		lrs = new LRSplineSurface(p+1, p+1, p+1, p+1);
	}

	// setup refinement parameters
	enum refinementStrategy strat;
	if(scheme == 0)
		strat = LR_FULLSPAN;
	else if(scheme == 1)
		strat = LR_MINSPAN;
	else if(scheme == 2)
		strat = LR_STRUCTURED_MESH;
	if(vol) {
		lrv->setRefMultiplicity(m);
		lrv->setRefStrat(strat);
	} else {
		lrs->setRefMultiplicity(m);
		lrs->setRefStrat(strat);
	}

	// for all iterations
	for(int n=0; n<N; n++) {
		vector<int> indices;

		if(scheme < 2) {
			if(vol) {
				lrv->generateIDs();
				lrv->getDiagonalElements(indices);
				lrv->refineElement(indices);
			} else {
				lrs->generateIDs();
				lrs->getDiagonalElements(indices);
				lrs->refineElement(indices);
			}
		} else if(scheme == 2) {
			if(vol) {
				// lrv->generateIDs();
				// lrv->getDiagonalBasisfunctions(indices);
				// lrv->refineBasisFunction(indices);
			} else {
				lrs->generateIDs();
				lrs->getDiagonalBasisfunctions(indices);
				lrs->refineBasisFunction(indices);
			}
		}

		// draw result files (with next refinement-step-diagonal shaded)
		if(dumpfile && !vol) {
			vector<int> diagElms, diagFuncs;
			lrs->getDiagonalElements(diagElms);
			lrs->getDiagonalBasisfunctions(diagFuncs);
			char filename[128];

			sprintf(filename, "index_%02d.eps", n+1);
			ofstream out;
			out.open(filename);
			if(scheme < 2)
				lrs->writePostscriptMesh(out, true, &diagElms);
			else
				lrs->writePostscriptMesh(out);
			out.close();

			sprintf(filename, "parameter_%02d.eps", n+1);
			out.open(filename);
			if(scheme < 2)
				lrs->writePostscriptElements(out, 2,2, true, &diagElms);
			else 
				lrs->writePostscriptElements(out, 2,2, true);
			out.close();

			sprintf(filename, "func_%02d.eps", n+1);
			out.open(filename);
			if(scheme == 2)
				lrs->writePostscriptFunctionSpace(out, &diagFuncs);
			else
				lrs->writePostscriptFunctionSpace(out);
			out.close();
			
			out.open("dof.m", ios::app);
			out << lrs->nBasisFunctions() << endl;
			out.close();

			out.open("elements.m", ios::app);
			out << lrs->nElements() << endl;
			out.close();
		}
	}
	
	vector<int> overloadedBasis;
	vector<int> overloadedElements;
	vector<int> multipleOverloadedElements;

	// harvest some statistics and display these results
	int nBasis      = 0;
	int nElements   = 0;
	int nMeshlines  = 0;
	if(vol) lrv->generateIDs();
	else    lrs->generateIDs();
	double avgBasisToElement = 0;
	double avgBasisToLine    = 0;
	int maxBasisToElement    = -1;
	int minBasisToElement    = 9999999;
	int nOverloadedElms      = 0;
	int nOverloadedBasis     = 0;
	if(vol) {
		nBasis     = lrv->nBasisFunctions();
		nElements  = lrv->nElements();
		nMeshlines = lrv->nMeshRectangles();
	} else {
		nBasis     = lrs->nBasisFunctions();
		nElements  = lrs->nElements();
		nMeshlines = lrs->nMeshlines();
	}
	HashSet<Basisfunction*> basis    = (vol) ? lrv->getAllBasisfunctions() : lrs->getAllBasisfunctions();
	vector<Element*>        elements = (vol) ? lrv->getAllElements()       : lrs->getAllElements();
	for(Basisfunction* b : basis) {
		int nE = b->nSupportedElements();
		maxBasisToElement = (maxBasisToElement > nE) ? maxBasisToElement : nE;
		minBasisToElement = (minBasisToElement < nE) ? minBasisToElement : nE;
		avgBasisToElement += nE;
		if(b->isOverloaded()) {
			nOverloadedBasis++;
			overloadedBasis.push_back(b->getId());
		}
	}
	avgBasisToElement /= nBasis;
	avgBasisToLine    /= nBasis;

	double avgElementToBasis       = 0;
	double avgSquareElementToBasis = 0;
	double hashCodePercentage      = ((double)basis.uniqueHashCodes())/nBasis;
	int maxElementToBasis          = -1;
	int minElementToBasis          = 9999999;
	for(Element *e : elements) {
		int nB = e->nBasisFunctions();
		maxElementToBasis       = (maxElementToBasis > nB) ? maxElementToBasis : nB;
		minElementToBasis       = (minElementToBasis < nB) ? minElementToBasis : nB;
		avgElementToBasis       += nB;
		avgSquareElementToBasis += nB*nB;
		if(e->isOverloaded()) {
			nOverloadedElms++;
			overloadedElements.push_back(e->getId());
		}
	}
	avgElementToBasis /= nElements;
	avgSquareElementToBasis /= nElements;
	
	cout << "Some statistics: " << endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Number of basisfunctions: " << nBasis           << endl;
	cout << "Number of elements      : " << nElements        << endl;
	cout << "Number of meshlines     : " << nMeshlines       << endl;
	cout << "-------------------------------------------------------------" << endl;
	cout << "Min number of Basisfuntion -> Element: " << minBasisToElement << endl;
	cout << "Max number of Basisfuntion -> Element: " << maxBasisToElement << endl;
	cout << "Avg number of Basisfuntion -> Element: " << avgBasisToElement << endl;
	cout << endl;
	cout << "Min number of        Element -> Basisfunction: " << minElementToBasis  << endl;
	cout << "Max number of        Element -> Basisfunction: " << maxElementToBasis  << endl;
	cout << "Avg number of        Element -> Basisfunction: " << avgElementToBasis  << endl;
	cout << "Avg square number of Element -> Basisfunction: " << avgSquareElementToBasis
	     << " (" << sqrt(avgSquareElementToBasis) << ")" << endl;
	cout << endl;
	cout << "Number of overloaded Basisfunctions : " << nOverloadedBasis  << endl;
	cout << "Number of overloaded Elements       : " << nOverloadedElms   << endl;
	cout << endl;
	cout << "Number of unique hashcodes          : " << basis.uniqueHashCodes() ;
	cout <<                                      " (" << hashCodePercentage*100 << " %)"  << endl;
	cout << "-------------------------------------------------------------" << endl;
#ifdef HAS_BOOST
	if(nBasis < 1300 && !vol) {
	cout << "Is linearly independent : " << ((lrs->isLinearIndepByMappingMatrix(false) )? "True":"False") << endl;
	cout << "-------------------------------------------------------------" << endl;
	}
#endif


	if(dumpfile) {
		cout << endl;
		cout << "Written";
		if(!vol) {
			ofstream meshfile;
			meshfile.open("mesh.eps");
			lrs->writePostscriptMesh(meshfile);
			meshfile.close();
			cout << " mesh to mesh.eps and";
		}
		
		ofstream lrfile;
		lrfile.open("diagonal.lr");
		if(vol) lrfile << *lrv << endl;
		else    lrfile << *lrs << endl;
		lrfile.close();
		
		cout << " diagonal.lr\n";
	}
}

