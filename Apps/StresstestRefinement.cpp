#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/MeshRectangle.h"
#include "LRSpline/HashSet.h"

using namespace LR;
using namespace std;

int main(int argc, char **argv) {
#ifdef TIME_LRSPLINE
	Profiler prof(argv[0]);
#endif

	/* UNIFORM refinement
	 * inserts global meshlines in a tensor product way, just stored within a 
	 * LRSplineSurface object
	 * 
	 * CORNER refinement
	 * splits the bottom left corner by introducing 3 new corner functions for
	 * each successively inserted cross. The knotlines are distributed uniformely
	 * throughout the lower left element
	 */
	enum refinement_scheme {UNIFORM, CORNER, DIAGONAL} refinement_scheme;

	// set default parameter values
	int goalBasisFunctions = 15000;
	int p1 = 3;
	int p2 = 3;
	int p3 = 3;
	int n1 = 7;
	int n2 = 7;
	int n3 = 7;
	int dim = 4;
	bool rat = false;
	bool vol = false;
	bool dumpFile = false;
	refinement_scheme = UNIFORM;

	string parameters(" parameters: \n" \
	                  "   -p1   <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	                  "   -p2   <n>  polynomial order in second parametric direction\n" \
	                  "   -p3   <n>  polynomial order in third parametric direction (volumes only)\n" \
	                  "   -p    <n>  polynomial order in all parametric directions\n" \
	                  "   -n1   <n>  number of basis functions in first parametric direction\n" \
	                  "   -n2   <n>  number of basis functions in second parametric direction\n" \
	                  "   -n3   <n>  number of basis functions in third parametric direction (volumes only)\n" \
	                  "   -n    <n>  number of basis functions in all parametric directions\n" \
	                  "   -dim  <n>  dimension of the controlpoints\n" \
	                  "   -goal <n>  number of basis functions before terminating program\n"\
	                  "   -unif      UNIFORM refinemen scheme\n"\
	                  "   -corner    CORNER refinemen scheme\n"\
	                  "   -diag      DIAGONAL refinemen scheme\n"\
	                  "   -vol       create a LRSplineVolume instead of Surface\n"\
	                  "   -dumpfile  writes an eps- and txt-file of the LR-mesh\n"\
	                  "   -help      display (this) help screen\n"\
	                  " <refine inputfile>\n"\
	                  "   inputfile describing meshline insertions\n");
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p3") == 0)
			p3 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p") == 0)
			p1 = p2 = p3 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n") == 0)
			n1 = n2 = n3 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n3") == 0)
			n3 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-dim") == 0)
			dim = atoi(argv[++i]);
		else if(strcmp(argv[i], "-goal") == 0)
			goalBasisFunctions = atoi(argv[++i]);
		else if(strcmp(argv[i], "-unif") == 0)
			refinement_scheme = UNIFORM;
		else if(strcmp(argv[i], "-corner") == 0)
			refinement_scheme = CORNER;
		else if(strcmp(argv[i], "-diag") == 0)
			refinement_scheme = DIAGONAL;
		else if(strcmp(argv[i], "-vol") == 0)
			vol = true;
		else if(strcmp(argv[i], "-dumpfile") == 0)
			dumpFile = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cerr << "usage: " << argv[0] << endl << parameters;
			exit(0);
		} else {
			cerr << "usage: " << argv[0] << endl << parameters;
			exit(1);
		}
	}

	// do some error testing on input
	if(n1 < p1) {
		cerr << "ERROR: n1 must be greater or equal to p1\n";
		exit(2);
	} else if(n2 < p2) {
		cerr << "ERROR: n2 must be greater or equal to p2\n";
		exit(2);
	} else if(n3 < p3) {
		cerr << "ERROR: n3 must be greater or equal to p3\n";
		exit(2);
	}

	// make a uniform integer knot vector
	double knot_u[n1+p1];
	double knot_v[n2+p2];
	double knot_w[n3+p3];
	int nCP = (vol) ? n1*n2*n3 : n1*n2;
	nCP    *= (dim+rat);
	double cp[nCP];
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;
	for(int i=0; i<p3+n3; i++)
		knot_w[i] = (i<p3) ? 0 : (i>n3) ? n3-p3+1 : i-p3+1;

	// create a list of random control points (all between 0.1 and 1.1)
	int k=0;
	for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
		cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0
		
	// make two identical surfaces
	LRSplineVolume  *lv;
	LRSplineSurface *lr;
	if(vol)
		lv = new LRSplineVolume (n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
	else                         
		lr = new LRSplineSurface(n1, n2,     p1, p2,     knot_u, knot_v,         cp, dim, rat);

	// perform the actual refinement
	int end1 = n1-p1+1;
	int end2 = n2-p2+1;
	int end3 = n3-p3+1;
	double h = 1.0;
	int nBasis      = 0;
	int nElements   = 0;
	int nMeshlines  = 0;
	int steppingDir = 0;
	double u = h/2.0;
	double v = h/2.0;
	double w = h/2.0;
	double unif_step_h = 1.0 / ((goalBasisFunctions - nBasis) / 3.0 + 1.0);
	int iter = 0;
	if(vol) {
		nBasis     = lv->nBasisFunctions();
		nElements  = lv->nElements();
		nMeshlines = lv->nMeshRectangles();
	} else {
		nBasis     = lr->nBasisFunctions();
		nElements  = lr->nElements();
		nMeshlines = lr->nMeshlines();
	}
	while(nBasis < goalBasisFunctions) {
		if(refinement_scheme == UNIFORM) {
			if(steppingDir == 0) {
				if(vol) lv->insert_line(new MeshRectangle(u,0,0,u,end2,end3));
				else    lr->insert_const_u_edge(u, 0, end2);
				u += h;
				if(u > end1) {
					steppingDir = 1;
					u = h/4.0;
				}
			} else if(steppingDir == 1) {
				if(vol) lv->insert_line(new MeshRectangle(0,v,0,end1,v,end3));
				else    lr->insert_const_v_edge(v, 0, end1);
				v += h;
				if(v > end2) {
					if(vol) {
						steppingDir = 2;
					} else {
						steppingDir = 0;
						h /= 2.0;
					}
					v = h/4.0;
				}
			} else if(steppingDir == 2) {
				lv->insert_line(new MeshRectangle(0,0,w,end1,end2,w));
				w += h;
				if(w > end3) {
					steppingDir = 0;
					w = h/4.0;
					h /= 2.0;
				}
			}
		} else if(refinement_scheme == CORNER) {
			if(vol) {
				lv->insert_line(new MeshRectangle(h-unif_step_h,0,0,  h-unif_step_h,h,h));
				lv->insert_line(new MeshRectangle(0,h-unif_step_h,0,  h,h-unif_step_h,h));
				lv->insert_line(new MeshRectangle(0,0,h-unif_step_h,  h,h,h-unif_step_h));
			} else {
				lr->insert_const_u_edge(h-unif_step_h, 0, h);
				lr->insert_const_v_edge(h-unif_step_h, 0, h);
			}
			h -= unif_step_h;
		} else if(refinement_scheme == DIAGONAL) {
			// guess this is hand tailored for degree p=3
			if(vol) {
				double startU = (iter-1<0)        ?    0 : (iter-1)*h;
				double stopU  = ((iter+2)*h>end1) ? end1 : (iter+2)*h;
				double startV = (iter-1<0)        ?    0 : (iter-1)*h;
				double stopV  = ((iter+2)*h>end2) ? end2 : (iter+2)*h;
				double startW = (iter-1<0)        ?    0 : (iter-1)*h;
				double stopW  = ((iter+2)*h>end3) ? end3 : (iter+2)*h;
				lv->insert_line(new MeshRectangle(  u   , startV, startW,   u  , stopV, stopW));
				lv->insert_line(new MeshRectangle(startU,   v   , startW, stopU,   v  , stopW));
				lv->insert_line(new MeshRectangle(startU, startV,   w   , stopU, stopV,   w  ));
			} else {
				lr->insert_const_u_edge(u, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end2) ? end2 : (iter+2)*h);
				lr->insert_const_v_edge(v, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end1) ? end1 : (iter+2)*h);
			}
			u += h;
			v += h;
			w += h;
			iter++;
			if( u>end1 ) {
				h /= 2.0;
				iter = 0;
				u = h/2.0;
				v = h/2.0;
				w = h/2.0;
			}
		}

		if(vol) {
			nBasis     = lv->nBasisFunctions();
			nElements  = lv->nElements();
			nMeshlines = lv->nMeshRectangles();
		} else {
			nBasis     = lr->nBasisFunctions();
			nElements  = lr->nElements();
			nMeshlines = lr->nMeshlines();
		}

		system("clear");
		cout << "LR type                  : " << ((vol)?"Volume":"Surface") << endl;
		cout << "Refinement scheme        : ";
		if(refinement_scheme == UNIFORM) cout << "UNIFORM\n";
		else if(refinement_scheme == CORNER) cout << "CORNER\n";
		else if(refinement_scheme == DIAGONAL) cout << "DIAGONAL\n";
		cout << "GOAL basis functions     : " << goalBasisFunctions << endl;
		cout << "=================================================" << endl;
		cout << endl;
		cout << "Number of basis functions: " << nBasis             << endl;
		cout << "Number of elements       : " << nElements          << endl;
		cout << "Number of meshlines      : " << nMeshlines         << endl;
	}

	// harvest some statistics and display these results
	double avgBasisToElement = 0;
	double avgBasisToLine    = 0;
	int maxBasisToElement    = -1;
	int minBasisToElement    = 9999999;

	HashSet<Basisfunction*> basis    = (vol) ? lv->getAllBasisfunctions() : lr->getAllBasisfunctions();
	vector<Element*>        elements = (vol) ? lv->getAllElements()       : lr->getAllElements();
	for(Basisfunction *b : basis) {
		int nE = b->nSupportedElements();
		maxBasisToElement = (maxBasisToElement > nE) ? maxBasisToElement : nE;
		minBasisToElement = (minBasisToElement < nE) ? minBasisToElement : nE;
		avgBasisToElement += nE;
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
	}
	avgElementToBasis       /= nElements;
	avgSquareElementToBasis /= nElements;
	
	cout << "Some statistics: " << endl;
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
	cout << "Number of unique hashcodes           : " << basis.uniqueHashCodes() ;
	cout <<                                      " (" << hashCodePercentage*100 << " %)"  << endl;

	
	if(dumpFile) {
		cout << endl;
		cout << "Written";
		if(!vol) {
			ofstream meshfile;
			meshfile.open("mesh.eps");
			lr->writePostscriptMesh(meshfile);
			meshfile.close();
			cout << " mesh to mesh.eps and";
		}
		
		ofstream lrfile;
		lrfile.open("stresstest.lr");
		if(vol) lrfile << *lv << endl;
		else    lrfile << *lr << endl;
		lrfile.close();
		
		cout << " stresstest.lr\n";
	}
}

