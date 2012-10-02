#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <GoTools/geometry/SplineSurface.h>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"

using namespace Go;
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
	int n1 = 7;
	int n2 = 7;
	int dim = 4;
	bool rat = false;
	bool dumpFile = false;
	refinement_scheme = UNIFORM;

	string parameters(" parameters: \n" \
	                  "   -p1   <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	                  "   -p2   <n>  polynomial order in second parametric direction\n" \
	                  "   -n1   <n>  number of basis functions in first parametric direction\n" \
	                  "   -n2   <n>  number of basis functions in second parametric direction\n" \
	                  "   -dim  <n>  dimension of the controlpoints\n" \
	                  "   -goal <n>  number of basis functions before terminating program\n"\
	                  "   -unif      UNIFORM refinemen scheme\n"\
	                  "   -corner    CORNER refinemen scheme\n"\
	                  "   -diag      DIAGONAL refinemen scheme\n"\
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
		else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
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
	}

	// make a uniform integer knot vector
	double knot_u[n1+p1];
	double knot_v[n2+p2];
	double cp[(dim+rat)*n1*n2];
	for(int i=0; i<p1+n1; i++)
		knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
	for(int i=0; i<p2+n2; i++)
		knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;

	// create a list of random control points (all between 0.1 and 1.1)
	int k=0;
	for(int j=0; j<n2; j++) 
		for(int i=0; i<n1; i++) 
			for(int d=0; d<dim+rat; d++)
				cp[k++] = ((i*2+j*3+d*5) % 13) / 13.0 + 0.1;  // rational weights also random and thus we need >0
		
	// make two identical surfaces
	SplineSurface   ss(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
	LRSplineSurface lr(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);

	// perform the actual refinement
	int end1 = n1-p1+1;
	int end2 = n2-p2+1;
	double h = 1.0;
	int nBasis = lr.nBasisFunctions();
	bool step_u = true;
	double u = h/2.0;
	double v = h/2.0;
	double unif_step_h = 1.0 / ((goalBasisFunctions - nBasis) / 3.0 + 1.0);
	int iter = 0;
	while(nBasis < goalBasisFunctions) {
		if(refinement_scheme == UNIFORM) {
			if(step_u) {
				lr.insert_const_u_edge(u, 0, end2);
				u += h;
				if(u > end1) {
					step_u = !step_u;
					u = h/4.0;
				}
			} else {
				lr.insert_const_v_edge(v, 0, end1);
				v += h;
				if(v > end2) {
					step_u = !step_u;
					v = h/4.0;
					h /= 2.0;
				}
			}
		} else if(refinement_scheme == CORNER) {
			lr.insert_const_u_edge(h-unif_step_h, 0, h);
			lr.insert_const_v_edge(h-unif_step_h, 0, h);
			h -= unif_step_h;
		} else if(refinement_scheme == DIAGONAL) {
			lr.insert_const_u_edge(u, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end2) ? end2 : (iter+2)*h);
			lr.insert_const_v_edge(v, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end1) ? end1 : (iter+2)*h);
			u += h;
			v += h;
			iter++;
			if( u>end1 ) {
				h /= 2.0;
				iter = 0;
				u = h/2.0;
				v = h/2.0;
			}
		}
		nBasis = lr.nBasisFunctions();
		int nElements = lr.nElements();
		int nMeshlines = lr.nMeshlines();

		system("clear");
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
	vector<Basisfunction*>::iterator bit;
	double avgBasisToElement = 0;
	double avgBasisToLine    = 0;
	int maxBasisToElement    = -1;
	int minBasisToElement    = 9999999;
	int maxBasisToLine       = -1;
	int minBasisToLine       = 9999999;
	for(bit=lr.basisBegin(); bit!=lr.basisEnd(); bit++) {
		int nE = (*bit)->nSupportedElements();
		int nL = (*bit)->nPartialLines();
		maxBasisToElement = (maxBasisToElement > nE) ? maxBasisToElement : nE;
		minBasisToElement = (minBasisToElement < nE) ? minBasisToElement : nE;
		avgBasisToElement += nE;
		maxBasisToLine = (maxBasisToLine > nL) ? maxBasisToLine : nL;
		minBasisToLine = (minBasisToLine < nL) ? minBasisToLine : nL;
		avgBasisToLine += nL;
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
	}
	avgElementToBasis /= lr.nElements();
	avgSquareElementToBasis /= lr.nElements();
	
	cout << "Some statistics: " << endl;
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



// currently I think we're going to avoid any verifications that this is right

#if 0
	// compare function values on edges, knots and in between the knots
	// as well as all derivatives (up to first derivatives)
	vector<Point> lr_pts(3), ss_pts(3);
	vector<bool> assert_function, assert_POF; // POF = partition of unity
	vector<double> par_u_values;
	vector<double> par_v_values;
	vector<Element*>::iterator el;
	for(el=lr.elementBegin(); el!=lr.elementEnd(); el++) {
		// evaluate one midpoint, 4 edge midpoints and 4 corners of the element
		double umin = (*el)->umin();
		double vmin = (*el)->vmin();
		double umax = (*el)->umax();
		double vmax = (*el)->vmax();
		int startI = ( umin == lr.startparam_u() ) ? 0 : 1;
		int startJ = ( vmin == lr.startparam_v() ) ? 0 : 1;
		for(int ii=startI; ii<3; ii++) {
			for(int jj=startJ; jj<3; jj++) {
				double u = umin + ii*(umax-umin)/2.0;
				double v = vmin + jj*(vmax-vmin)/2.0;
				lr.point(lr_pts, u,v, 1);
				ss.point(ss_pts, u,v, 1, false, false);

				bool correct = true;
				for(int i=0; i<3; i++)  
					for(int d=0; d<dim; d++)
						if( fabs((lr_pts[i][d]-ss_pts[i][d])/ss_pts[i][d]) > TOL ) // relative error
							correct = false;
				assert_function.push_back(correct);
				par_u_values.push_back(u);
				par_v_values.push_back(v);

				// test for partition of unity
				BasisDerivsSf lr_basis;
				lr.computeBasis(u,v, lr_basis);
				double sum        = 0;
				double sum_diff_u = 0;
				double sum_diff_v = 0;
				for(uint i=0; i<lr_basis.basisValues.size(); i++) {
					sum        += lr_basis.basisValues[i] ;
					sum_diff_u += lr_basis.basisDerivs_u[i];
					sum_diff_v += lr_basis.basisDerivs_v[i];
				}
				cout << "   LR(" << u << ", " << v << ") = " << lr_pts[0] << endl;
				cout << "   SS(" << u << ", " << v << ") = " << ss_pts[0] << endl;
				cout << "dx LR(" << u << ", " << v << ") = " << lr_pts[1] << endl;
				cout << "dx SS(" << u << ", " << v << ") = " << ss_pts[1] << endl;
				cout << "dy LR(" << u << ", " << v << ") = " << lr_pts[2] << endl;
				cout << "dy SS(" << u << ", " << v << ") = " << ss_pts[2] << endl;
				cout << "sum (" << u << ", " << v << ") minus one = " << sum-1.0 << endl;
				cout << "sum diff u = " << sum_diff_u << endl;
				cout << "sum diff v = " << sum_diff_v << endl;
				
				correct = !(fabs(sum-1.0) > TOL || fabs(sum_diff_u) > TOL || fabs(sum_diff_v) > TOL);
				assert_POF.push_back(correct);
			}
		}
	}
	bool oneFail = false;
	bool linearIndep ;
	bool doActualLinTest ;
	if(lr.nBasisFunctions() > max_n_linear_depence_testing) {
		doActualLinTest = false;
		linearIndep = true;
	} else {
		doActualLinTest = true;
		linearIndep = lr.isLinearIndepByMappingMatrix(true);
	}
	cout << endl;
	cout << "   =====================    RESULT SUMMARY   ====================     \n\n";
	cout << "_____u____________v_____________FUNCTION______PARTITON OF UNITY____\n";
	for(uint i=0; i<assert_function.size(); i++) {
		printf("%10.4g   %10.4g           ", par_u_values[i], par_v_values[i]);
		if(assert_function[i])
			cout << "OK              ";
		else
			cout << "FAIL            ";
		if(assert_POF[i])
			cout << "OK\n";
		else
			cout << "FAIL\n";
		if(!assert_function[i])
			oneFail = true;
		if(!assert_POF[i])
			oneFail = true;
	}
	cout << "----------------------------------------\n";
	cout << "  Linear independent :     " << ((doActualLinTest) ? ((linearIndep) ? "OK" : "FAIL") : "(NOT TESTED)") << endl;
	cout << "----------------------------------------\n";
	if(oneFail || !linearIndep)
		cout << "    test FAILED\n";
	else
		cout << "    all assertions passed\n";
	cout << "========================================\n";
	cout << "\n\n\n";
	cout << "Key LR-spline information:\n";
	cout << "  number of basis functions: " << lr.nBasisFunctions() << endl;
	cout << "  number of mesh lines     : " << lr.nMeshlines() << endl;
	cout << "  number of elements       : " << lr.nElements() << endl;

#endif
	
	
	if(dumpFile) {
		ofstream meshfile;
		meshfile.open("mesh.eps");
		lr.writePostscriptMesh(meshfile);
		meshfile.close();
		
		ofstream lrfile;
		lrfile.open("lrspline.txt");
		lrfile << lr << endl;
		lrfile.close();
		
		cout << endl;
		cout << "Written mesh to mesh.eps and lrspline.txt\n";
	}
}

