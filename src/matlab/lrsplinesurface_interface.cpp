#include <cmath>
#include <iostream>
#include <sstream>
#include <LRSpline/LRSplineSurface.h>
#include <mex.h>
#include <matrix.h>
#include "class_handle.hpp"

using namespace std;
using namespace LR;

#define DOUBLE_TOL 1e-12

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	char cmd[64];
	if(nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string, less than 64 characters long");


	/*********************************************************
	 *******      Constructors and destructors             ***
	 ********************************************************/
	
	if(!strcmp(cmd, "new")) {
		// Check parameters
		if (nlhs != 1)
			mexErrMsgTxt("New: One output expected.");
		
		// fetch input arguments
		double *n = mxGetPr(prhs[1]);
		double *p = mxGetPr(prhs[2]);
		if(nrhs > 4) {
			double *knotU = mxGetPr(prhs[3]);
			double *knotV = mxGetPr(prhs[4]);
			if(nrhs == 6) {
				double *cp = mxGetPr(prhs[5]);
				int    dim = mxGetM(prhs[5]);

		// Return a handle to a new C++ instance
				plhs[0] = convertPtr2Mat<LRSplineSurface>(new LRSplineSurface(floor(n[0]), floor(n[1]), floor(p[0])+1, floor(p[1])+1, knotU, knotV, cp, dim, false));
			} else {
				plhs[0] = convertPtr2Mat<LRSplineSurface>(new LRSplineSurface(floor(n[0]), floor(n[1]), floor(p[0])+1, floor(p[1])+1, knotU, knotV));
			}
		} else {
			plhs[0] = convertPtr2Mat<LRSplineSurface>(new LRSplineSurface(floor(n[0]), floor(n[1]), floor(p[0])+1, floor(p[1])+1));
		}
		return;
	}

	// Check there is a second input, which should be the class instance handle
	if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
	
	// Delete
	if (!strcmp("delete", cmd)) {
		// Destroy the C++ object
		destroyObject<LRSplineSurface>(prhs[1]);
		// Warn if other commands were ignored
		if (nlhs != 0 || nrhs != 2)
			mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
		return;
	}


	// Get the class instance pointer from the second input
	LRSplineSurface *lr = convertMat2Ptr<LRSplineSurface>(prhs[1]);


	// Copy
	if (!strcmp("copy", cmd)) {
		// Warn if other commands were ignored
		if (nlhs != 1 || nrhs != 2)
			mexWarnMsgTxt("Copy: Unexpected arguments ignored.");
		plhs[0] = convertPtr2Mat<LRSplineSurface>(lr->copy());
		return;
	}
	

	/*********************************************************
	 *******      Call the various class methods           ***
	 ********************************************************/


	// Print    
	if (!strcmp("print", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 2)
			mexErrMsgTxt("Print: Unexpected arguments.");
		// write the object into a stringstream and rewrite this using mexPrintf
		stringstream ss;
		lr->write(ss);
		mexPrintf(ss.str().c_str());
		return;
	}


	// Get Element Containing
	if (!strcmp("get_element_containing", cmd)) {
		// Check parameters
		if (nlhs < 1 || nrhs < 3)
			mexErrMsgTxt("compute_basis: Unexpected arguments.");
		double *u = mxGetPr(prhs[2]);

		int id = lr->getElementContaining(u[0], u[1]); // always search for element and return sparse result
		plhs[0] = mxCreateDoubleScalar(id+1);
		return;
	}


	// Compuate Basisfunctions
	if (!strcmp("compute_basis", cmd)) {
		// Check parameters
		if (nlhs < 1 || nrhs < 3)
			mexErrMsgTxt("compute_basis: Unexpected arguments.");
		int derivs = 0;
		if(nrhs > 3)
			derivs = floor(mxGetScalar(prhs[3]));
		double *u = mxGetPr(prhs[2]);

		int id = lr->getElementContaining(u[0], u[1]); // always search for element and return sparse result

		vector<vector<double> > vResult;
		lr->computeBasis(u[0], u[1], vResult, derivs, id);
		plhs[0] = mxCreateDoubleMatrix(vResult[0].size(),vResult.size(), mxREAL);
		double *dResult = mxGetPr(plhs[0]);
		for(size_t i=0; i<vResult.size(); i++)
			copy(vResult[i].begin(), vResult[i].end(), dResult+i*vResult[0].size());
		
		return;
	}


	// Get Bezier Extraction Matrix
	if (!strcmp("get_bezier_extraction", cmd)) {
		// Check parameters
		if (nlhs < 1 || nrhs < 3)
			mexErrMsgTxt("compute_basis: Unexpected arguments.");

		int iEl          = floor(mxGetScalar(prhs[2])) - 1;
		Element *element = lr->getElement(iEl);

		vector<double> vResult;
		lr->getBezierExtraction(iEl, vResult);
		plhs[0] = mxCreateDoubleMatrix(element->nBasisFunctions(),lr->order(0)*lr->order(1), mxREAL);
		double *dResult = mxGetPr(plhs[0]);
		copy(vResult.begin(), vResult.end(), dResult);
		
		return;
	}


	// insert line
	if (!strcmp("insert_line", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 4)
			mexErrMsgTxt("insert_line: Unexpected arguments.");

		// rewrap results into vector array
		double *start = mxGetPr(prhs[2]);
		double *stop  = mxGetPr(prhs[3]);
		for(int i=2; i<4; i++) {
			int m = mxGetM(prhs[i]);
			int n = mxGetN(prhs[i]);
			if(m*n != 2)
				mexErrMsgTxt("insert_line: Invalid arguments");
		}
		int mult = 1;
		if(nrhs > 3)
			mult = floor(mxGetScalar(prhs[4]));

		if( fabs(start[0] - stop[0]) < DOUBLE_TOL)
			lr->insert_const_u_edge(start[0], start[1], stop[1], mult);
		else if( fabs(start[1] - stop[1]) < DOUBLE_TOL)
			lr->insert_const_v_edge(start[1], start[0], stop[0], mult);
		else
			mexErrMsgTxt("insert_line: Invalid arguments");

		return;
	}



	// Refine basisfunctions
	if (!strcmp("refine_basis", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 4)
			mexErrMsgTxt("Refine_basis: Unexpected arguments.");

		// rewrap results into vector array
		double *b = mxGetPr(prhs[2]);
		int m = mxGetM(prhs[2]);
		int n = mxGetN(prhs[2]);
		int mult = floor(mxGetScalar(prhs[3]));
		vector<int> basis(m*n);
		for(size_t i=0; i<n*m; i++)
			basis[i] = floor(b[i]-1);
		
		lr->setRefMultiplicity(mult);
		lr->refineBasisFunction(basis);
		return;
	}


	// Refine elements
	if (!strcmp("refine_elements", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 4)
			mexErrMsgTxt("Refine_elements: Unexpected arguments.");

		double *el = mxGetPr(prhs[2]);
		int m    = mxGetM(prhs[2]);
		int n    = mxGetN(prhs[2]);
		int mult = floor(mxGetScalar(prhs[3]));
		vector<int> elements(m*n);
		for(size_t i=0; i<n*m; i++)
			elements[i] = floor(el[i]-1);
		
		lr->setRefMultiplicity(mult);
		lr->refineElement(elements);
		return;
	}


	// Refine basisfunctions
	if (!strcmp("raise_order", cmd)) {
		// Check parameters
		if (nlhs < 1 || nrhs < 4)
			mexErrMsgTxt("Raise_order: Unexpected arguments.");

		int dp = floor(mxGetScalar(prhs[2]));
		int dq = floor(mxGetScalar(prhs[3]));
		
		plhs[0] = convertPtr2Mat<LRSplineSurface>(lr->getRaiseOrderSpace(dp, dq));
		return;
	}


	// Refine basisfunctions
	if (!strcmp("set_control_points", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 3)
			mexErrMsgTxt("Set_control_points: Unexpected arguments.");

		double *data = mxGetPr(prhs[2]);
		int     m    = mxGetM(prhs[2]);
		int     n    = mxGetN(prhs[2]);

		if(m != lr->dimension() || n != lr->nBasisFunctions())
			mexErrMsgTxt("Set_control_points: invalid controlpoint size ");

		vector<double> cps(m*n);
		copy(data, data+m*n, cps.begin());
		if(! lr->setControlPoints(cps) )
			mexErrMsgTxt("Set_control_points: C++ class failed. Unknown error");
		return;
	}


	// Set Continutiy
	if (!strcmp("set_continuity", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 3)
			mexErrMsgTxt("set_continuity: Unexpected arguments.");

		double *c = mxGetPr(prhs[2]);
		int     m    = mxGetM(prhs[2]);
		int     n    = mxGetN(prhs[2]);
		if(m*n != 2)
			mexErrMsgTxt("set_continuity: invalid continuity size ");

		lr->setGlobalContinuity(c[0], c[1]);
		return;
	}


	// Point
	if (!strcmp("point", cmd)) {
		// Check parameters
		if (nlhs < 1 || nrhs < 3)
			mexErrMsgTxt("Point: Unexpected arguments.");

		double *u = mxGetPr(prhs[2]);
		vector<double> vResult;
		lr->point(vResult, u[0], u[1]);
		plhs[0] = mxCreateDoubleMatrix(vResult.size(),1, mxREAL);
		double *dResult  = mxGetPr(plhs[0]);
		copy(vResult.begin(), vResult.end(), dResult);

		return;
	}


	// update primitives
	if (!strcmp("get_primitives", cmd)) {
		// Check parameters
		if (nlhs < 0 || nrhs < 2)
			mexErrMsgTxt("Get_primitives: Unexpected arguments.");
		int p1 = lr->order(0);
		int p2 = lr->order(1);
		lr->generateIDs();
		// write basis functions
		if(nlhs > 0) {
			int n   = lr->nBasisFunctions();
			int dim = lr->dimension();
			plhs[0]  = mxCreateDoubleMatrix(n, p1+p2+2, mxREAL);
			double *knots = mxGetPr(plhs[0]);
			double *cp;
			double *w;
			if(nlhs>1) {
				plhs[1]  = mxCreateDoubleMatrix(dim, n, mxREAL);
				cp       = mxGetPr(plhs[1]);
			}
			if(nlhs>2) {
				plhs[2]  = mxCreateDoubleMatrix(1, n, mxREAL);
				w        = mxGetPr(plhs[2]);
			}
			HashSet_iterator<Basisfunction*> it;
			int i=0;
			for(it=lr->basisBegin(); it!=lr->basisEnd(); ++it, ++i) {
				for(int j=0; j<p1+1; j++)
					knots[j*n+i] = (**it)[0][j]; // knot vector in u-direction
				for(int j=0; j<p2+1; j++)
					knots[(j+p1+1)*n+i] = (**it)[1][j]; // knot vector in v-direction
				if(nlhs>1)
					for(int j=0; j<dim; j++)
						cp[i*dim+j] = (**it).cp(j); // control points
				if(nlhs>2)
					w[i] = (**it).w(); // weights
			}
		}
		if(nlhs > 3) {
			int n = lr->nMeshlines();
			plhs[3] = mxCreateDoubleMatrix(n, 5, mxREAL);
			double *meshline = mxGetPr(plhs[3]);
			vector<Meshline*>::iterator it;
			int i=0;
			for(it=lr->meshlineBegin(); it!=lr->meshlineEnd(); ++it, ++i) {
				if((**it).span_u_line_) {
					meshline[i+0*n] = (**it).start_;
					meshline[i+1*n] = (**it).const_par_;
					meshline[i+2*n] = (**it).stop_;
					meshline[i+3*n] = (**it).const_par_;
				} else {
					meshline[i+0*n] = (**it).const_par_;
					meshline[i+1*n] = (**it).start_;
					meshline[i+2*n] = (**it).const_par_;
					meshline[i+3*n] = (**it).stop_;
				}
				meshline[i+4*n] = (**it).multiplicity_;
			}
		}
		if(nlhs > 4) {
			int n = lr->nElements();
			plhs[4] = mxCreateDoubleMatrix(n, 4, mxREAL);
			double *element = mxGetPr(plhs[4]);
			double *support = mxGetPr(plhs[4]);
			if(nlhs > 5) {
				plhs[5] = mxCreateCellMatrix(n, 1);
			}
			vector<Element*>::iterator it;
			int i=0;
			for(it=lr->elementBegin(); it!=lr->elementEnd(); ++it, ++i) {
				element[i+0*n] = (**it).getParmin(0);
				element[i+1*n] = (**it).getParmin(1);
				element[i+2*n] = (**it).getParmax(0);
				element[i+3*n] = (**it).getParmax(1);
				if(nlhs > 5) {
					int m = (**it).nBasisFunctions();
					mxArray *arraySupport = mxCreateDoubleMatrix(1, m, mxREAL);
					double  *elSupport = mxGetPr(arraySupport);
					HashSet_iterator<Basisfunction*> bit;
					int j=0;
					for(bit=(**it).supportBegin(); bit!=(**it).supportEnd(); ++bit, ++j)
						elSupport[j] = (**bit).getId() + 1;
					mxSetCell(plhs[5], i, arraySupport);
				}
			}
		}
		if(nlhs > 6) {
			int n = lr->nElements();
			plhs[6] = mxCreateDoubleMatrix(2, 1, mxREAL);
			double *p = mxGetPr(plhs[6]);
			p[0] = lr->order(0)-1;
			p[1] = lr->order(1)-1;
		}
		return;
	}



	/*********************************************************
	 *****      Error handling: no commands recognized   *****
	 ********************************************************/

	// Got here, so command not recognized
	mexErrMsgTxt("Command not recognized.");

}
