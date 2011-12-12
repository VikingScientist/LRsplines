
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
	int p      = 3;
	int N      = 6;
	int m      = 1;
	int scheme = 0;
	string parameters(" parameters: \n" \
	                  "   -p      <n> polynomial DEGREE (order-1) of the basis\n" \
	                  "   -n      <n> number of iterations\n" \
	                  "   -m      <n> knot multiplicity\n" \
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
	ofstream clear_out;
	clear_out.open("dof.m");
	clear_out.close();
	clear_out.open("elements.m");
	clear_out.close();



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
	
	for(int n=0; n<N; n++) {
		lr.generateIDs();
		vector<int> indices;
		if(scheme < 2) {
			vector<Element*>::iterator iel;
			for(iel=lr.elementBegin(); iel<lr.elementEnd(); iel++) {
				if((**iel).umin() == (**iel).vmin())
					indices.push_back((**iel).getId() );
			}
		} else if(scheme == 2) {
			vector<Basisfunction*>::iterator ifun;
			for(ifun=lr.basisBegin(); ifun<lr.basisEnd(); ifun++) {
				bool diag = true;
				for(int i=0; i<lr.order_u()+1; i++)
					diag = diag && (**ifun).knot_u_[i] == (**ifun).knot_v_[i];
				if(diag)
					indices.push_back( (**ifun).getId() );
			}
		}
		if(scheme == 0)
			lr.refineElement(indices, m, false, false);
		else if(scheme == 1)
			lr.refineElement(indices, m, true,  false);
		else if(scheme == 2)
			lr.refineBasisFunctions(indices, m);

		char filename[128];

		sprintf(filename, "index_%02d.eps", n+1);
		ofstream out;
		out.open(filename);
		lr.writePostscriptMesh(out, true, (scheme < 2));
		out.close();

		sprintf(filename, "parameter_%02d.eps", n+1);
		out.open(filename);
		lr.writePostscriptElements(out, true, (scheme < 2));
		out.close();

		sprintf(filename, "func_%02d.eps", n+1);
		out.open(filename);
		lr.writePostscriptFunctionSpace(out, (scheme == 2));
		out.close();
		
		out.open("dof.m", ios::app);
		out << lr.nBasisFunctions() << endl;
		out.close();

		out.open("elements.m", ios::app);
		out << lr.nElements() << endl;
		out.close();
	}

}

