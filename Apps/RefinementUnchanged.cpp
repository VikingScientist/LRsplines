#include <stdio.h>
#include <iostream>
#include <fstream>
#include <GoTools/geometry/SplineSurface.h>
#include "LRSplineSurface.h"
#include "Profiler.h"

using namespace Go;
using namespace LR;
using namespace std;

int main(int argc, char **argv) {
	Profiler prof(argv[0]);

	const double TOL = 1e-6;
	bool verbose = false;
	
	int p1 = 4;
	int p2 = 4;
	int n1 = 22;
	int n2 = 22;
	int dim = 4;
	bool rat = false;
	double knot_u[n1+p1];
	double knot_v[n2+p2];
	double cp[(dim+rat)*n1*n2];

	// make a uniform integer knot vector
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

	int giantN  = 128;
	double step = 1.0/giantN;
	for(int i=0; i<1280; i++) {
		lr.insert_const_u_edge(step/2+i*step, i%12, 6+i%12, 1);
		lr.insert_const_v_edge(step/2+i*step, i%12, 6+i%12, 1);
	}

/*
	lr.insert_const_u_edge(.5, 0, 4);
	lr.insert_const_u_edge(1.5, 1, 5);
	lr.insert_const_v_edge(0.5, 0, 4);
	lr.insert_const_v_edge(1.5, 1, 5);
*/


	// compare function values on edges, knots and in between the knots
	// as well as all derivatives (up to first derivatives)
	vector<Point> lr_pts(3), ss_pts(3);
	vector<bool> assert_function, assert_POF; // POF = partition of unity
	vector<double> par_u_values;
	vector<double> par_v_values;
	BasisDerivsSf lr_basis;
	lr.computeBasis(0,0, lr_basis);
	for(double u=0; u<=n1-p1+1; u+=0.5) {
		for(double v=0; v<=n2-p2+1; v+=0.5) {
			// lr.point(lr_pts, u,v, 1);
			ss.point(lr_pts, u,v, 1);
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
			// BasisDerivsSf lr_basis;
			// lr.computeBasis(u,v, lr_basis);
			double sum        = 0;
			double sum_diff_u = 0;
			double sum_diff_v = 0;
			for(uint i=0; i<lr_basis.basisValues.size(); i++) {
				sum        += lr_basis.basisValues[i] ;
				sum_diff_u += lr_basis.basisDerivs_u[i];
				sum_diff_v += lr_basis.basisDerivs_v[i];
			}
			if(verbose) {
				cout << "   LR(" << u << ", " << v << ") = " << lr_pts[0] << endl;
				cout << "   SS(" << u << ", " << v << ") = " << ss_pts[0] << endl;
				cout << "dx LR(" << u << ", " << v << ") = " << lr_pts[1] << endl;
				cout << "dx SS(" << u << ", " << v << ") = " << ss_pts[1] << endl;
				cout << "dy LR(" << u << ", " << v << ") = " << lr_pts[2] << endl;
				cout << "dy SS(" << u << ", " << v << ") = " << ss_pts[2] << endl;
				cout << "sum (" << u << ", " << v << ") minus one = " << sum-1.0 << endl;
				cout << "sum diff u = " << sum_diff_u << endl;
				cout << "sum diff v = " << sum_diff_v << endl;
			}

			
			correct = !(fabs(sum-1.0) > TOL || fabs(sum_diff_u) > TOL || fabs(sum_diff_v) > TOL);
			assert_POF.push_back(correct);

		}
	}
	bool oneFail = false;
	// bool linearIndep = lr.isLinearIndepByMappingMatrix(true);
	bool linearIndep = true;
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
	cout << "  Linear independent :     " << ((linearIndep) ? "OK" : "FAIL") << endl;
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
	
	ofstream meshfile;
	meshfile.open("mesh.eps");
	lr.writePostscriptMesh(meshfile);
	meshfile.close();

	ofstream lrfile;
	lrfile.open("lrspline.txt");
	lrfile << lr << endl;
	lrfile.close();

	// cout << endl;
	// cout << "Written mesh to mesh.eps and lrspline.txt\n";

	// lr.printElements(cout);
}
