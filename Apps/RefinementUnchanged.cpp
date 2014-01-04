#include <cmath>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#ifdef HAS_GOTOOLS
	#include <GoTools/geometry/SplineSurface.h>
	#include <GoTools/trivariate/SplineVolume.h>
	#include <GoTools/geometry/ObjectHeader.h>
#endif
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/MeshRectangle.h"

typedef unsigned int uint;

using namespace LR;
using namespace std;

ostream& operator<<(ostream& out, const vector<double> vec) {
	for(double d : vec)
		out << d << " ";
	return out;
}

int main(int argc, char **argv) {
#ifdef TIME_LRSPLINE
	Profiler prof(argv[0]);
#endif

	// set default parameter values
	const double TOL = 1e-6;
	const double max_n_linear_depence_testing = 1000;
	int p1 = 3;
	int p2 = 2;
	int p3 = 4;
	int n1 = 6;
	int n2 = 5;
	int n3 = 8;
	double *knot_u      = NULL;
	double *knot_v      = NULL;
	double *knot_w      = NULL;
	int dim             = 4;
	int nDiagonals      = -1;
	bool rat            = false;
	bool dumpFile       = false;
	bool vol            = false;
	char *lrInitMesh    = NULL;
	char *inputFileName = NULL;
	stringstream parameters;
	parameters << " parameters: \n" \
	              "   -p1    <n>  polynomial ORDER (degree+1) in first parametric direction\n" \
	              "   -p2    <n>  polynomial order in second parametric direction\n" \
	              "   -p3    <n>  polynomial order in third parametric direction (only trivariate volumes)\n" \
	              "   -n1    <n>  number of basis functions in first parametric direction\n" \
	              "   -n2    <n>  number of basis functions in second parametric direction\n" \
	              "   -n3    <n>  number of basis functions in third parametric direction (only trivariate volumes)\n" \
	              "   -knot1 <n>  space-seperated list of the first knot vector (must specify n1 and p1 first)\n"\
	              "   -knot2 <n>  space-seperated list of the second knot vector (must specify n2 and p2 first)\n"\
	              "   -knot3 <n>  space-seperated list of the third knot vector (must specify n3 and p3 first)\n"\
	              "   -dim   <n>  dimension of the controlpoints\n" \
	              "   -diag  <n>  override inputfile and run diagonal testcase\n"\
	              "   -in:   <s>  make the LRSplineSurface <s> the initial mesh\n"\
	              "   -dumpfile   writes an eps- and txt-file of the LR-mesh (bivariate surfaces only)\n"\
	              "   -help       display (this) help screen\n";
	parameters << " default values\n";
	parameters << "   -p   = { " << p1 << ", " << p2 << ", " << p3 << " }\n";
	parameters << "   -n   = { " << n1 << ", " << n2 << ", " << n3 << " }\n";
	parameters << "   -vol = { " << ((vol)?"true":"false") << " }\n";
	parameters << " <refine inputfile>\n"\
	              "   inputfile describing meshline insertions.\n"\
	              "   FORMAT:\n"\
	              "     <numb. inserted lines>\n"\
	              "     <is_const_u> <const_par> <start> <stop> <mult>                    (Surface only)\n"\
	              "     <constParDir> <constPar> <start1> <stop1> <start2> <stop2> <mult> (Volume only)\n";
	
	// read input
	for(int i=1; i<argc; i++) {
		if(strcmp(argv[i], "-p1") == 0)
			p1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p2") == 0)
			p2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-p3") == 0) {
			p3  = atoi(argv[++i]);
			vol = true;
		} else if(strcmp(argv[i], "-n1") == 0)
			n1 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n2") == 0)
			n2 = atoi(argv[++i]);
		else if(strcmp(argv[i], "-n3") == 0) {
			n3  = atoi(argv[++i]);
			vol = true;
		}  else if(strcmp(argv[i], "-dim") == 0)
			dim = atoi(argv[++i]);
		else if(strncmp(argv[i], "-in:", 4) == 0)
			lrInitMesh = argv[i]+4;
		else if(strcmp(argv[i], "-diag") == 0)
			nDiagonals = atoi(argv[++i]);
		else if(strcmp(argv[i], "-vol") == 0)
			vol = true;
		else if(strcmp(argv[i], "-dumpfile") == 0)
			dumpFile = true;
		else if(strcmp(argv[i], "-help") == 0) {
			cout << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters.str();
			exit(0);
		} else if(strcmp(argv[i], "-knot1") == 0) {
			knot_u = new double[n1+p1];
			for(int j=0; j<n1+p1; j++)
				knot_u[j] = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-knot2") == 0) {
			knot_v = new double[n2+p2];
			for(int j=0; j<n2+p2; j++)
				knot_v[j] = atoi(argv[++i]);
		} else if(strcmp(argv[i], "-knot3") == 0) {
			knot_w = new double[n3+p3];
			for(int j=0; j<n3+p3; j++)
				knot_w[j] = atoi(argv[++i]);
		} else {
			if(inputFileName != NULL) {
				cerr << "usage: " << argv[0] << endl << parameters.str();
				exit(1);
			} else {
				inputFileName = argv[i];
			}
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
	} else if(nDiagonals==-1 && inputFileName == NULL) {
		cerr << "ERROR: Specify input file name\n";
		cerr << "usage: " << argv[0] << "[parameters] <refine inputfile>" << endl << parameters.str();
		exit(3);
	}

#ifdef HAS_GOTOOLS
	Go::SplineSurface   *ss;
	Go::SplineVolume    *sv;
#else 
	LRSplineSurface *ss;
	LRSplineVolume  *sv;
#endif
	LRSplineSurface *lrs;
	LRSplineVolume  *lrv;

	if(lrInitMesh == NULL) {

		// make a uniform integer knot vector
		if(knot_u == NULL) {
			knot_u = new double[n1+p1];
			for(int i=0; i<p1+n1; i++)
				knot_u[i] = (i<p1) ? 0 : (i>n1) ? n1-p1+1 : i-p1+1;
		}
		if(knot_v == NULL) {
			knot_v = new double[n2+p2];
			for(int i=0; i<p2+n2; i++)
				knot_v[i] = (i<p2) ? 0 : (i>n2) ? n2-p2+1 : i-p2+1;
		}
		if(knot_w == NULL) {
			knot_w = new double[n3+p3];
			for(int i=0; i<p3+n3; i++)
				knot_w[i] = (i<p3) ? 0 : (i>n3) ? n3-p3+1 : i-p3+1;
		}

		// create a list of random control points (all between 0.1 and 1.1)
		int nCP = (vol) ? n1*n2*n3 : n1*n2;
		nCP    *= (dim+rat);
		double cp[nCP];
		int k=0;
		for(int i=0; i<nCP; i++) // 839 as a generator over Z_853 gives a period of 425. Should suffice
			cp[k++] = (i*839 % 853) / 853.0 + 0.1;  // rational weights also random and thus we need >0
			
		// make two spline objects (using GoTools for reference if available)
		if(vol) {
#ifdef HAS_GOTOOLS
			sv  = new Go::SplineVolume(n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
#else 
			sv  = new LRSplineVolume  (n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
#endif
			lrv = new LRSplineVolume  (n1, n2, n3, p1, p2, p3, knot_u, knot_v, knot_w, cp, dim, rat);
		} else {
#ifdef HAS_GOTOOLS
			ss  = new Go::SplineSurface(n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
#else 
			ss  = new LRSplineSurface  (n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
#endif
			lrs = new LRSplineSurface  (n1, n2, p1, p2, knot_u, knot_v, cp, dim, rat);
		}
	} else {
#ifdef HAS_GOTOOLS
		ifstream inputfile;
		inputfile.open(lrInitMesh);
		if(!inputfile.is_open()) {
			cerr << "Error: could not open file " << lrInitMesh << endl;
			exit(3);
		}
		Go::ObjectHeader head;
		ss = new Go::SplineSurface();
		sv = new Go::SplineVolume();
		inputfile >> head;
		if(head.classType() == Go::Class_SplineVolume) {
			vol = true;
			inputfile >> *sv;
			lrv = new LRSplineVolume(sv);
			dim = sv->dimension();
			rat = sv->rational();
			n1  = sv->numCoefs(0);
			n2  = sv->numCoefs(1);
			n3  = sv->numCoefs(2);
			p1  = sv->order(0);
			p2  = sv->order(1);
			p3  = sv->order(2);
		} else if(head.classType() == Go::Class_SplineSurface) {
			vol = false;
			inputfile >> *ss;
			lrs = new LRSplineSurface(ss);
			dim = ss->dimension();
			rat = ss->rational();
			n1  = ss->numCoefs_u();
			n2  = ss->numCoefs_v();
			p1  = ss->order_u();
			p2  = ss->order_v();
		}
#else
		cerr << "Error: could not open file " << lrInitMesh << endl;
		exit(4);
#endif
	}

	vector<bool>   is_const_u;  // for surfaces
	vector<int>    constParDir; // for volumes
	vector<double> constPar;
	vector<double> startPar1;
	vector<double> endPar1;
	vector<double> startPar2;
	vector<double> endPar2;
	vector<int> multiplicity;
	if(nDiagonals==-1) {
		// read input-file 
		ifstream inputFile;
		inputFile.open(inputFileName);
		if( inputFile.is_open() ) {
			int n;
			inputFile >> n;
			int f;
			double a,b,c,d,e;
			int m;
			for(int i=0; i<n; i++) {
				inputFile >> f;
				inputFile >> a;
				inputFile >> b;
				inputFile >> c;
				if(vol) {
					inputFile >> d;
					inputFile >> e;
				}
				inputFile >> m;

				is_const_u.push_back(f); // these are used for surfaces
				constPar.push_back(a);
				startPar1.push_back(b);
				endPar1.push_back(c);
				multiplicity.push_back(m);

				constParDir.push_back(f); // these are used for volumes
				startPar2.push_back(d);
				endPar2.push_back(e);
			}
			inputFile.close();
		} else {
			cerr <<"ERROR: could not open file " << inputFileName << endl;
			exit(3);
		}
	}


	if(nDiagonals==-1) {
		if(vol) {
			for(uint i=0; i<constParDir.size(); i++) {
				int c  = constParDir[i]; // constant dir
				int v1 = (c+1)%3;        // first variable dir
				int v2 = (c+2)%3;        // second variable dir
				double min[3];
				double max[3];
				min[c]  = constPar[i];
				min[v1] = startPar1[i];
				min[v2] = startPar2[i];
				max[c]  = constPar[i];
				max[v1] = endPar1[i];
				max[v2] = endPar2[i];

				lrv->insert_line(new MeshRectangle(min[0], min[1], min[2], max[0], max[1], max[2]));
			}
		} else {
			for(uint i=0; i<is_const_u.size(); i++) {
				if(is_const_u[i])
					lrs->insert_const_u_edge(constPar[i], startPar1[i], endPar1[i], multiplicity[i]);
				else
					lrs->insert_const_v_edge(constPar[i], startPar1[i], endPar1[i], multiplicity[i]);
			}
		}
	} else {
		for(int diagRuns=0; diagRuns<nDiagonals; diagRuns++) {
			vector<Element*>::iterator it;
			vector<int> diagonalElements;
			int i=0;
			for(it=lrs->elementBegin(); it<lrs->elementEnd(); it++, i++)
				if((**it).umin() == (**it).vmin())
					diagonalElements.push_back(i);
			lrs->refineElement(diagonalElements);
		}
	}
	

	// compare function values on edges, knots and in between the knots
	// as well as all derivatives (up to first derivatives)
	bool oneFail         = false;
	bool linearIndep     = true;
	bool doActualLinTest ;
	vector<bool> assert_function, assert_POF; // POF = partition of unity
	if(!vol) {
		/***************************************************************************
		******                SURFACE TESTING                                  *****
		****************************************************************************/
#ifdef HAS_GOTOOLS
		vector<Go::Point> lr_pts(3), ss_pts(3);
#else
		vector<vector<double> > lr_pts(3), ss_pts(3);
#endif
		vector<double> par_u_values;
		vector<double> par_v_values;
		vector<Element*>::iterator el;
		int el_i=0;
		for(el=lrs->elementBegin(); el!=lrs->elementEnd(); el++, el_i++) {
			// evaluate one midpoint, 4 edge midpoints and 4 corners of the element
			double umin = (*el)->umin();
			double vmin = (*el)->vmin();
			double umax = (*el)->umax();
			double vmax = (*el)->vmax();
			int startI = ( umin == lrs->startparam(0) ) ? 0 : 1;
			int startJ = ( vmin == lrs->startparam(1) ) ? 0 : 1;
			for(int ii=startI; ii<3; ii++) {
				for(int jj=startJ; jj<3; jj++) {
					double u = umin + ii*(umax-umin)/2.0;
					double v = vmin + jj*(vmax-vmin)/2.0;
					lrs->point(lr_pts, u,v, 1);
					ss->point(ss_pts, u,v, 1);

					bool correct = true;
					for(int i=0; i<3; i++)  
						for(int d=0; d<dim; d++)
							if( fabs(lr_pts[i][d]-ss_pts[i][d]) > TOL ) // relative error
								correct = false;
					assert_function.push_back(correct);
					par_u_values.push_back(u);
					par_v_values.push_back(v);

					// test for partition of unity
					vector<vector<double> > lr_basis;
					lrs->computeBasis(u,v, lr_basis, 1);
					double sum        = 0;
					double sum_diff_u = 0;
					double sum_diff_v = 0;
					for(uint i=0; i<lr_basis.size(); i++) {
						sum        += lr_basis[i][0];
						sum_diff_u += lr_basis[i][1];
						sum_diff_v += lr_basis[i][2];
					}
					cout << "   LR(" << u << ", " << v << ") = " << lr_pts[0] << endl;
					cout << "   SS(" << u << ", " << v << ") = " << ss_pts[0] << endl;
					cout << "dx LR(" << u << ", " << v << ") = " << lr_pts[1] << endl;
					cout << "dx SS(" << u << ", " << v << ") = " << ss_pts[1] << endl;
					cout << "dy LR(" << u << ", " << v << ") = " << lr_pts[2] << endl;
					cout << "dy SS(" << u << ", " << v << ") = " << ss_pts[2] << endl;
					cout << "sum ("  << u << ", " << v << ") minus one = " << sum-1.0 << endl;
					cout << "sum diff u = " << sum_diff_u << endl;
					cout << "sum diff v = " << sum_diff_v << endl;
					
					correct = !(fabs(sum-1.0) > TOL || fabs(sum_diff_u) > TOL || fabs(sum_diff_v) > TOL);
					assert_POF.push_back(correct);
				}
			}
		}
		if(lrs->nBasisFunctions() > max_n_linear_depence_testing) {
			doActualLinTest = false;
			linearIndep = true;
		} else {
			doActualLinTest = true;
			linearIndep = lrs->isLinearIndepByMappingMatrix(true);
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
	} else {
		/***************************************************************************
		******                VOLUME  TESTING                                  *****
		****************************************************************************/
#ifdef HAS_GOTOOLS
		vector<Go::Point> lr_pts(4), ss_pts(4);
#else
		vector<vector<double> > lr_pts(4), ss_pts(4);
#endif
		vector<double> par_u_values;
		vector<double> par_v_values;
		vector<double> par_w_values;
		vector<Element*>::iterator el;
		int el_i=0;
		for(el=lrv->elementBegin(); el!=lrv->elementEnd(); el++, el_i++) {
			// evaluate one midpoint, 4 edge midpoints and 4 corners of the element
			double umin = (*el)->getParmin(0);
			double vmin = (*el)->getParmin(1);
			double wmin = (*el)->getParmin(2);
			double umax = (*el)->getParmax(0);
			double vmax = (*el)->getParmax(1);
			double wmax = (*el)->getParmax(2);
			int startI = ( umin == lrv->startparam(0) ) ? 0 : 1;
			int startJ = ( vmin == lrv->startparam(1) ) ? 0 : 1;
			int startK = ( wmin == lrv->startparam(2) ) ? 0 : 1;
			for(int ii=startI; ii<3; ii++) {
				for(int jj=startJ; jj<3; jj++) {
					for(int kk=startK; kk<3; kk++) {
						double u = umin + ii*(umax-umin)/2.0;
						double v = vmin + jj*(vmax-vmin)/2.0;
						double w = wmin + kk*(wmax-wmin)/2.0;
						lrv->point(lr_pts, u,v,w, 1);
						sv->point( ss_pts, u,v,w, 1);

						bool correct = true;
						for(int i=0; i<3; i++)  
							for(int d=0; d<dim; d++)
								if( fabs(lr_pts[i][d]-ss_pts[i][d]) > TOL ) // relative error
									correct = false;
						assert_function.push_back(correct);
						par_u_values.push_back(u);
						par_v_values.push_back(v);
						par_w_values.push_back(w);

						// test for partition of unity
						vector<vector<double> > lr_basis;
						lrv->computeBasis(u,v,w, lr_basis, 1);
						double sum        = 0;
						double sum_diff_u = 0;
						double sum_diff_v = 0;
						double sum_diff_w = 0;
						for(uint i=0; i<lr_basis.size(); i++) {
							sum        += lr_basis[i][0];
							sum_diff_u += lr_basis[i][1];
							sum_diff_v += lr_basis[i][2];
							sum_diff_w += lr_basis[i][3];
						}
						cout << "   LR(" << u << ", " << v << ", " << w << ") = " << lr_pts[0] << endl;
						cout << "   SS(" << u << ", " << v << ", " << w << ") = " << ss_pts[0] << endl;
						cout << "dx LR(" << u << ", " << v << ", " << w << ") = " << lr_pts[1] << endl;
						cout << "dx SS(" << u << ", " << v << ", " << w << ") = " << ss_pts[1] << endl;
						cout << "dy LR(" << u << ", " << v << ", " << w << ") = " << lr_pts[2] << endl;
						cout << "dy SS(" << u << ", " << v << ", " << w << ") = " << ss_pts[2] << endl;
						cout << "dz LR(" << u << ", " << v << ", " << w << ") = " << lr_pts[3] << endl;
						cout << "dz SS(" << u << ", " << v << ", " << w << ") = " << ss_pts[3] << endl;
						cout << "sum ("  << u << ", " << v << ", " << w << ") minus one = " << sum-1.0 << endl;
						cout << "sum diff u = " << sum_diff_u << endl;
						cout << "sum diff v = " << sum_diff_v << endl;
						cout << "sum diff w = " << sum_diff_w << endl;
						
						correct = !(fabs(sum-1.0)    > TOL ||
						            fabs(sum_diff_u) > TOL ||
						            fabs(sum_diff_v) > TOL ||
						            fabs(sum_diff_w) > TOL  );
						assert_POF.push_back(correct);
					}
				}
			}
		}
		if(lrv->nBasisFunctions() > max_n_linear_depence_testing) {
			doActualLinTest = false;
			linearIndep = true;
		} else {
			// messed up because the mapping matrix is not operational yet. Will fix in due time
			doActualLinTest = true;
			doActualLinTest = false; 
			// linearIndep = lrv->isLinearIndepByMappingMatrix(true);
		}
		cout << endl;
		cout << "   =====================    RESULT SUMMARY   ====================     \n\n";
		cout << "_____u____________v____________w_____________FUNCTION______PARTITON OF UNITY____\n";
		for(uint i=0; i<assert_function.size(); i++) {
			printf("%10.4g   %10.4g   %10.4g           ", par_u_values[i], par_v_values[i], par_w_values[i]);
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


	int nBasis;
	int nMeshlines;
	int nElements;
	if(vol) {
		nBasis     = lrv->nBasisFunctions();
		nMeshlines = lrv->nMeshRectangles();
		nElements  = lrv->nElements();
	} else {
		nBasis     = lrs->nBasisFunctions();
		nMeshlines = lrs->nMeshlines();
		nElements  = lrs->nElements();
	}


	cout << "Key LR-spline information:\n";
	cout << "  number of basis functions: " << nBasis     << endl;
	cout << "  number of mesh lines     : " << nMeshlines << endl;
	cout << "  number of elements       : " << nElements  << endl;
	
	
	if(dumpFile) {
		cout << endl;
		cout << "Written ";

		if(!vol) {
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
			cout << "mesh.eps, functions.eps, domain.eps, controlmesh.eps and ";
		}
		cout << "RefUnchagned.lr\n";
	
		ofstream lrfile;
		lrfile.open("RefUnchagned.lr");
		if(vol) lrfile << *lrv << endl;
		else    lrfile << *lrs << endl;
		lrfile.close();
	}

}
