
#include <GoTools/trivariate/SplineVolume.h>
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"
#include <iostream>
#include <fstream>
#include <cstring>

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {

	// set default parameter values
	int p         = 3;
	int N         = 3;
	int m         = 1;
	int scheme    = 0;  
	double R      = 0.5;   // radius of refinement layer
	double eps    = 0.005; // layer width parameter
	string parameters(" parameters: \n" \
	                  "   -p      <n> polynomial DEGREE (order-1) of the basis\n" \
	                  "   -n      <n> number of iterations\n" \
	                  "   -m      <n> knot multiplicity\n" \
	                  "   -scheme <n> refinement scheme (0=FULLSPAN, 1=MINSPAN)\n" \
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
	if(scheme > 1)
		scheme = 1;
	if(scheme < 0)
		scheme = 0;

	// setup refinement parameters
	enum refinementStrategy strat;
	if(scheme == 0)
		strat = LR_SAFE;
	else if(scheme == 1)
		strat = LR_MINSPAN;

	// make a uniform integer knot vector
	double knot_u[] = {0,0,1,1};
	double knot_v[] = {0,0,1,1};
	double knot_w[] = {0,0,1,1};
	double cp3[]    = {0,0,0,
	                   1,0,0,
	                   0,1,0,
	                   1,1,0,
	                   0,0,1,
	                   1,0,1,
	                   0,1,1,
	                   1,1,1};


	// create initial geometry and an empty value thing
	SplineVolume   sv(2, 2, 2,  2, 2, 2, knot_u, knot_v, knot_w, cp3, 3, false); // WHY does this not compile!?
	sv.raiseOrder(p-1, p-1, p-1);
	LRSplineVolume *lrGeom   = new LRSplineVolume(&sv);
	LRSplineVolume *lrValues = lrGeom->copy();

	// do a variational diminishing approximation on the LRSpline object
	HashSet_iterator<Basisfunction*> geom = lrGeom->basisBegin();
	lrValues->rebuildDimension(1);
	for(Basisfunction *b : lrValues->getAllBasisfunctions()) {
		double x = (**geom).cp(0);
		double y = (**geom).cp(1);
		double z = (**geom).cp(2);
		double r = sqrt(x*x + y*y + z*z);
		*b->cp() = atan((r-0.5)/eps);
		++geom;
	}

	// make room for all iterations... geometries and field values
	vector<LRSplineVolume*> geometries;
	vector<LRSplineVolume*> values;
	geometries.push_back(lrGeom);
	values.push_back(lrValues);


	// for all iterations
	for(int n=0; n<N; n++) {
		std::cout << "Starting iteration " << (n+1) << " / " << N <<  std::endl;
		lrGeom = geometries.back()->copy();
		lrGeom->setRefMultiplicity(m);
		lrGeom->setRefStrat(strat);
		lrGeom->generateIDs();

		// find all elements which intersect the solution layer
		vector<int> indices;
		int i=0;
		for(Element *el : lrGeom->getAllElements()) {
			double x0 = el->getParmin(0);
			double y0 = el->getParmin(1);
			double z0 = el->getParmin(2);
			double x1 = el->getParmax(0);
			double y1 = el->getParmax(1);
			double z1 = el->getParmax(2);
			if(x0*x0 + y0*y0 + z0*z0 < R*R &&
			   x1*x1 + y1*y1 + z1*z1 > R*R) {
				indices.push_back(el->getId());
			}
			i++;
		}

		// perform actual refinement
		lrGeom->refineElement(indices);

		// make a VDA approximation to the solution
		lrValues = lrGeom->copy();
		geom = lrGeom->basisBegin();
		lrValues->rebuildDimension(1);
		for(Basisfunction *b : lrValues->getAllBasisfunctions()) {
			double x = (**geom).cp(0);
			double y = (**geom).cp(1);
			double z = (**geom).cp(2);
			double r = sqrt(x*x + y*y + z*z);
			*b->cp() = atan((r-0.5)/eps)  /  b->w();
			++geom;
		}

		// add to final storage
		values.push_back(lrValues);
		geometries.push_back(lrGeom);
	}
	
	
	// write results to file
	ofstream lrfile;
	lrfile.open("waterfall_geom.lr");
	lrfile << *geometries.back() << endl;
	lrfile.close();
	lrfile.open("waterfall_values.lr");
	lrfile << *values.back() << endl;
	lrfile.close();
	cout << "Written waterfall_geom.lr and waterfall_values.lr" << endl;
}

