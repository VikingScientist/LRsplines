#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <boost/rational.hpp>
#include <cfloat>
#include <cmath>

typedef unsigned int uint;

//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;

namespace LR {

#define DOUBLE_TOL 1e-14
#define MY_STUPID_FABS(x) (((x)>0)?(x):-(x))


LRSplineSurface::LRSplineSurface() {
	rational_ = false;
	dim_      = 0;
	order_u_  = 0;
	order_v_  = 0;
	start_u_  = 0;
	start_v_  = 0;
	end_u_    = 0;
	end_v_    = 0;
	meshline_ = std::vector<Meshline*>(0);
	element_  = std::vector<Element*>(0);
	maxTjoints_       = -1;
	doCloseGaps_      = true;
	maxAspectRatio_   = 2.0;
	doAspectRatioFix_ = false;
	refStrat_         = LR_SAFE;
	refKnotlineMult_  = 1;
	symmetry_         = 1;
	element_red       = 0.5;
	element_green     = 0.5;
	element_blue      = 0.5;
	basis_red         = 1.0;
	basis_green       = 0.83;
	basis_blue        = 0.0;
	selected_basis_red    = 1.0;
	selected_basis_green  = 0.2;
	selected_basis_blue   = 0.05;
}

LRSplineSurface::LRSplineSurface(Go::SplineSurface *surf) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	rational_ = surf->rational();
	dim_      = surf->dimension();
	order_u_  = surf->order_u();
	order_v_  = surf->order_v();
	start_u_  = surf->startparam_u();
	start_v_  = surf->startparam_v();
	end_u_    = surf->endparam_u();
	end_v_    = surf->endparam_v();
	maxTjoints_       = -1;
	doCloseGaps_      = true;
	maxAspectRatio_   = 2.0;
	doAspectRatioFix_ = false;
	refStrat_         = LR_SAFE;
	refKnotlineMult_  = 1;
	symmetry_         = 1;
	element_red       = 0.5;
	element_green     = 0.5;
	element_blue      = 0.5;
	basis_red         = 1.0;
	basis_green       = 0.83;
	basis_blue        = 0.0;
	selected_basis_red    = 1.0;
	selected_basis_green  = 0.2;
	selected_basis_blue   = 0.05;


	int n1 = surf->numCoefs_u();
	int n2 = surf->numCoefs_v();
	std::vector<double>::iterator coef  = surf->coefs_begin();
	std::vector<double>::const_iterator knot_u = surf->basis(0).begin();
	std::vector<double>::const_iterator knot_v = surf->basis(1).begin();
	for(int j=0; j<n2; j++)
		for(int i=0; i<n1; i++)
			basis_.insert(new Basisfunction(knot_u+i, knot_v+j, coef+(j*n1+i)*(dim_+rational_), dim_, order_u_, order_v_));
	int unique_u=0;
	int unique_v=0;
	for(int i=0; i<n1+order_u_; i++) {// const u, spanning v
		int mult = 1;
		while(i+1<n1+order_u_ && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2], mult) );
	}
	for(int i=0; i<n2+order_v_; i++) {// const v, spanning u
		int mult = 1;
		while(i+1<n2+order_v_ && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n1], mult) );
	}
	for(int j=0; j<unique_v-1; j++) {
		for(int i=0; i<unique_u-1; i++) {
			double umin = meshline_[i]->const_par_;
			double vmin = meshline_[unique_u + j]->const_par_;
			double umax = meshline_[i+1]->const_par_;
			double vmax = meshline_[unique_u + j+1]->const_par_;
			element_.push_back(new Element(umin, vmin, umax, vmax));
		}
	}
	for(Basisfunction* b : basis_)
		updateSupport(b);
}


LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	rational_ = rational;
	dim_      = dim;
	order_u_  = order_u;
	order_v_  = order_v;
	start_u_  = knot_u[0];
	start_v_  = knot_v[0];
	end_u_    = knot_u[n1];
	end_v_    = knot_v[n2];
	maxTjoints_       = -1;
	doCloseGaps_      = true;
	maxAspectRatio_   = 2.0;
	doAspectRatioFix_ = false;
	refStrat_         = LR_SAFE;
	refKnotlineMult_  = 1;
	symmetry_         = 1;
	element_red       = 0.5;
	element_green     = 0.5;
	element_blue      = 0.5;
	basis_red         = 1.0;
	basis_green       = 0.83;
	basis_blue        = 0.0;
	selected_basis_red    = 1.0;
	selected_basis_green  = 0.2;
	selected_basis_blue   = 0.05;

	for(int j=0; j<n2; j++)
		for(int i=0; i<n1; i++)
			basis_.insert(new Basisfunction(knot_u+i, knot_v+j, coef+(j*n1+i)*(dim+rational), dim, order_u, order_v));
	int unique_u=0;
	int unique_v=0;
	for(int i=0; i<n1+order_u; i++) {// const u, spanning v
		int mult = 1;
		while(i<n1+order_u && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2], mult) );
	}
	for(int i=0; i<n2+order_v; i++) {// const v, spanning u
		int mult = 1;
		while(i<n2+order_v && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n1], mult) );
	}
	for(int j=0; j<unique_v-1; j++) {
		for(int i=0; i<unique_u-1; i++) {
			double umin = meshline_[i]->const_par_;
			double vmin = meshline_[unique_u + j]->const_par_;
			double umax = meshline_[i+1]->const_par_;
			double vmax = meshline_[unique_u + j+1]->const_par_;
			element_.push_back(new Element(umin, vmin, umax, vmax));
		}
	}

	for(Basisfunction* b : basis_)
		updateSupport(b);
}

LRSplineSurface::~LRSplineSurface() {
	for(Basisfunction* b : basis_)
		delete b;
	for(uint i=0; i<meshline_.size(); i++)
		delete meshline_[i];
	for(uint i=0; i<element_.size(); i++)
		delete element_[i];
}


LRSplineSurface* LRSplineSurface::copy() const {
	generateIDs();

	std::vector<Basisfunction*> basisVector;

	// flat list to make it quicker to update pointers from Basisfunction to Element and back again
	LRSplineSurface *returnvalue = new LR::LRSplineSurface();
	
	for(Basisfunction* b : basis_) {
		Basisfunction *newB = b->copy();
		returnvalue -> basis_.insert(newB);
		basisVector.push_back(newB);
	}
	
	for(Element *e : element_) {
		Element* newEl = e->copy();
		returnvalue -> element_.push_back(newEl);
		newEl->updateBasisPointers(basisVector);
	}
	
	for(Meshline *m : meshline_)
		returnvalue -> meshline_.push_back(m->copy());
	
	returnvalue->rational_         = this->rational_;
	returnvalue->dim_              = this->dim_;
	returnvalue->order_u_          = this->order_u_;
	returnvalue->order_v_          = this->order_v_;
	returnvalue->start_u_          = this->start_u_;
	returnvalue->start_v_          = this->start_v_;
	returnvalue->end_u_            = this->end_u_;
	returnvalue->end_v_            = this->end_v_;
	returnvalue->maxTjoints_       = this->maxTjoints_;
	returnvalue->doCloseGaps_      = this->doCloseGaps_;
	returnvalue->doAspectRatioFix_ = this->doAspectRatioFix_;
	returnvalue->maxAspectRatio_   = this->maxAspectRatio_;
	
	return returnvalue;
}

void LRSplineSurface::point(Go::Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const {
	Go::Point cp;
	double basis_ev;

	pt.resize(dim_);
	pt.setValue(0.0);
	
	for(Basisfunction* b : element_[iEl]->support()) {
		b->getControlPoint(cp);
		
		basis_ev = b->evaluate(u,v, u_from_right, v_from_right);
		pt += basis_ev*cp;
	}
	
}

void LRSplineSurface::point(Go::Point &pt, double u, double v, int iEl) const {
	Go::Point cp;
	double basis_ev;

	pt.resize(dim_);
	pt.setValue(0.0);
	if(iEl == -1)
		iEl = getElementContaining(u,v);

	for(Basisfunction* b : element_[iEl]->support()) {
		b->getControlPoint(cp);
		
		basis_ev = b->evaluate(u,v, u!=end_u_, v!=end_v_);
		pt += basis_ev*cp;
	}
	
}

void LRSplineSurface::point(std::vector<Go::Point> &pts, double u, double v, int derivs, int iEl) const {
#ifdef TIME_LRSPLINE
	PROFILE("Point()");
#endif
	Go::Point cp;
	std::vector<double> basis_ev;

	// clear and resize output array (optimization may consider this an outside task)
	pts.resize((derivs+1)*(derivs+2)/2);
	for(uint i=0; i<pts.size(); i++) {
		pts[i].resize(dim_);
		pts[i].setValue(0.0);
	}

	if(iEl == -1)
		iEl = getElementContaining(u,v);
	for(Basisfunction* b : element_[iEl]->support() ) {
		b->getControlPoint(cp);
		b->evaluate(basis_ev, u,v, derivs, u!=end_u_, v!=end_v_);
		for(uint j=0; j<pts.size(); j++)
			pts[j] += basis_ev[j]*cp;
	}
}

void LRSplineSurface::computeBasis (double param_u, double param_v, Go::BasisDerivsSf2 & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	std::vector<double> values;
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.prepareDerivs(param_u, param_v, 0, -1, nPts);

	int i=0;
	//element_[i]->write(std::cout);
	
	for(it=itStart; it!=itStop; ++it, ++i) {
		(*it)->evaluate(values, param_u, param_v, 2, param_u!=end_u_, param_v!=end_v_);
	
		result.basisValues[i]    = values[0];
		result.basisDerivs_u[i]  = values[1];
		result.basisDerivs_v[i]  = values[2];
		result.basisDerivs_uu[i] = values[3];
		result.basisDerivs_uv[i] = values[4];
		result.basisDerivs_vv[i] = values[5];
	}
}

void LRSplineSurface::computeBasis (double param_u, double param_v, Go::BasisDerivsSf & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	std::vector<double> values;
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.prepareDerivs(param_u, param_v, 0, -1, nPts);
	
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i) {
		(*it)->evaluate(values, param_u, param_v, 1, param_u!=end_u_, param_v!=end_v_);
		
		result.basisValues[i]   = values[0];
		result.basisDerivs_u[i] = values[1];
		result.basisDerivs_v[i] = values[2];
	}
}


void LRSplineSurface::computeBasis(double param_u, double param_v, Go::BasisPtsSf & result, int iEl ) const {
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();

	result.preparePts(param_u, param_v, 0, -1, nPts);
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i)
		result.basisValues[i] = (*it)->evaluate(param_u, param_v, param_u!=end_u_, param_v!=end_v_);
}

void LRSplineSurface::computeBasis (double param_u,
                                    double param_v,
                                    std::vector<std::vector<double> >& result,
                                    int derivs,
                                    int iEl ) const
{
	result.clear();
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	
	result.resize(nPts);
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i)
	    (*it)->evaluate(result[i], param_u, param_v, derivs, param_u!=end_u_, param_v!=end_v_);

}

int LRSplineSurface::getElementContaining(double u, double v) const {
	for(uint i=0; i<element_.size(); i++)
		if(element_[i]->umin() <= u && element_[i]->vmin() <= v) 
			if((u < element_[i]->umax() || (u == end_u_ && u <= element_[i]->umax())) && 
			   (v < element_[i]->vmax() || (v == end_v_ && v <= element_[i]->vmax())))
				return i;
		
	return -1;
}

void LRSplineSurface::getMinspanLines(int iEl, std::vector<Meshline*>& lines) {
	Element *e = element_[iEl];
	double umin = e->umin();
	double umax = e->umax();
	double vmin = e->vmin();
	double vmax = e->vmax();
	double min_du = DBL_MAX;
	double min_dv = DBL_MAX;
	int    best_startI = order_u_+2;
	int    best_stopI  = order_u_+2;
	int    best_startJ = order_v_+2;
	int    best_stopJ  = order_v_+2;
	bool   only_insert_span_u_line = (vmax-vmin) >= maxAspectRatio_*(umax-umin);
	bool   only_insert_span_v_line = (umax-umin) >= maxAspectRatio_*(vmax-vmin);
	// loop over all supported B-splines and choose the minimum one
	for(Basisfunction* b : e->support()) {
		double lowu  = b->getParmin(0);
		double highu = b->getParmax(0);
		double lowv  = b->getParmin(1);
		double highv = b->getParmax(1);
		double du = highu - lowu;
		double dv = highv - lowv;
		int startI=0;
		int stopI=0;
		int startJ=0;
		int stopJ=0;
		while((*b)[0][startI] <= e->umin())
			startI++;
		while((*b)[0][stopI]  <  e->umax())
			stopI++;
		while((*b)[1][startJ] <= e->vmin())
			startJ++;
		while((*b)[1][stopJ]  <  e->vmax())
			stopJ++;

		// min_du is defined as the minimum TOTAL knot span (of an entire basis function)
		bool fixU = false;
		int delta_startI = abs(startI - (order_u_+1)/2);
		int delta_stopI  = abs(stopI  - (order_u_+1)/2);
		if(  du <  min_du )
			fixU = true;
		if( du == min_du && delta_startI <= best_startI && delta_stopI  <= best_stopI )
			fixU = true;
		if(fixU) {
			umin = lowu;
			umax = highu;
			min_du = umax-umin;
			best_startI = delta_startI;
			best_stopI  = delta_stopI;
		}
		bool fixV = false;
		int delta_startJ = abs(startJ - (order_v_+1)/2);
		int delta_stopJ  = abs(stopJ  - (order_v_+1)/2);
		if(  dv <  min_dv )
			fixV = true;
		if( dv == min_dv && delta_startJ <= best_startJ && delta_stopJ  <= best_stopJ )
			fixV = true;
		if(fixV) {
			vmin = lowv;
			vmax = highv;
			min_dv = vmax-vmin;
			best_startJ = delta_startJ;
			best_stopJ  = delta_stopJ;
		}
	}

	if(!only_insert_span_v_line) 
		lines.push_back(new Meshline(true, (e->vmin() + e->vmax())/2.0, umin, umax, 1));
		
	if(!only_insert_span_u_line) 
		lines.push_back(new Meshline(false, (e->umin() + e->umax())/2.0, vmin, vmax, 1));

}

void LRSplineSurface::getFullspanLines(int iEl, std::vector<Meshline*>& lines) {
	Element *e = element_[iEl];
	double umin = e->umin();
	double umax = e->umax();
	double vmin = e->vmin();
	double vmax = e->vmax();
	bool   only_insert_span_u_line = (vmax-vmin) >= maxAspectRatio_*(umax-umin);
	bool   only_insert_span_v_line = (umax-umin) >= maxAspectRatio_*(vmax-vmin);
	// loop over all supported B-splines and make sure that everyone is covered by meshline
	for(Basisfunction *b : e->support()) {
		umin = (umin > (*b).getParmin(0)) ? (*b).getParmin(0) : umin;
		umax = (umax < (*b).getParmax(0)) ? (*b).getParmax(0) : umax;
		vmin = (vmin > (*b).getParmin(1)) ? (*b).getParmin(1) : vmin;
		vmax = (vmax < (*b).getParmax(1)) ? (*b).getParmax(1) : vmax;
	}
	if(!only_insert_span_v_line) 
		lines.push_back(new Meshline(true, (e->vmin() + e->vmax())/2.0, umin, umax, 1));
		
	if(!only_insert_span_u_line) 
		lines.push_back(new Meshline(false, (e->umin() + e->umax())/2.0, vmin, vmax, 1));
}

void LRSplineSurface::getStructMeshLines(int iBasis, std::vector<Meshline*>& lines) {
	Basisfunction *b = getBasisfunction(iBasis);
	double umin = b->getParmin(0);
	double umax = b->getParmax(0);
	double vmin = b->getParmin(1);
	double vmax = b->getParmax(1);

	// find the largest knotspan in this function
	double max_du = 0;
	double max_dv = 0;
	for(int j=0; j<order_u_; j++) {
		double du = (*b)[0][j+1]-(*b)[0][j];
		bool isZeroSpan =  MY_STUPID_FABS(du) < DOUBLE_TOL ;
		max_du = (isZeroSpan || max_du>du) ? max_du : du;
	}
	for(int j=0; j<order_v_; j++) {
		double dv = (*b)[1][j+1]-(*b)[1][j];
		bool isZeroSpan =  MY_STUPID_FABS(dv) < DOUBLE_TOL ;
		max_dv = (isZeroSpan || max_dv>dv) ? max_dv : dv;
	}

	// to keep as "square" basis function as possible, only insert
	// into the largest knot spans
	for(int j=0; j<order_u_; j++) {
		double du = (*b)[0][j+1]-(*b)[0][j];
		if( MY_STUPID_FABS(du-max_du) < DOUBLE_TOL )
			lines.push_back(new Meshline(false, ((*b)[0][j] + (*b)[0][j+1])/2.0, vmin, vmax,1));
	}
	for(int j=0; j<order_v_; j++) {
		double dv = (*b)[1][j+1]-(*b)[1][j];
		if( MY_STUPID_FABS(dv-max_dv) < DOUBLE_TOL )
			lines.push_back(new Meshline(true, ((*b)[1][j] + (*b)[1][j+1])/2.0, umin, umax,1));
	}
}

void LRSplineSurface::refineBasisFunction(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineBasisFunction(tmp);
}

void LRSplineSurface::refineBasisFunction(const std::vector<int> &indices) {
	std::vector<Meshline*> newLines;

	/* first retrieve all meshlines needed */
	for(uint i=0; i<indices.size(); i++)
		getStructMeshLines(indices[i],newLines);

	/* Do the actual refinement */
	for(uint i=0; i<newLines.size(); i++) {
		Meshline *m = newLines[i];
		insert_line(!m->is_spanning_u(), m->const_par_, m->start_, m->stop_, refKnotlineMult_);
	}

	/* do a posteriori fixes to ensure a proper mesh */
	aPosterioriFixes();

	/* exit cleanly be deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++) 
		delete newLines[i];
}

void LRSplineSurface::refineElement(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineElement(tmp);
}

void LRSplineSurface::refineElement(const std::vector<int> &indices) {
	std::vector<Meshline*> newLines;

	/* first retrieve all meshlines needed */
	for(uint i=0; i<indices.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanLines(indices[i],newLines);
		else
			getFullspanLines(indices[i],newLines);
	}

	/* Do the actual refinement */
	for(uint i=0; i<newLines.size(); i++) {
		Meshline *m = newLines[i];
		insert_line(!m->is_spanning_u(), m->const_par_, m->start_, m->stop_, refKnotlineMult_);
	}

	/* do a posteriori fixes to ensure a proper mesh */
	aPosterioriFixes();

	/* exit cleanly be deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++) 
		delete newLines[i];
}

void LRSplineSurface::refineByDimensionIncrease(const std::vector<double> &errPerElement, double beta) {
	Element       *e;
	/* accumulate the error & index - vector */
	std::vector<IndexDouble> errors;
	if(refStrat_ == LR_ISOTROPIC_FUNC) { // error per-function
		int i=0;
		for(Basisfunction *b : basis_) {
			errors.push_back(IndexDouble(0.0, i));
			for(int j=0; j<b->nSupportedElements(); j++) {
				e = *(b->supportedElementBegin() + j);
				errors[i].first += errPerElement[e->getId()];
			}
			i++;
		}
	} else {
		for(uint i=0; i<element_.size(); i++) 
			errors.push_back(IndexDouble(errPerElement[i], i));
	}

	/* sort errors */
	std::sort(errors.begin(), errors.end(), std::greater<IndexDouble>());

	/* first retrieve all possible meshlines needed */
	std::vector<std::vector<Meshline*> > newLines(errors.size(), std::vector<Meshline*>(0));
	for(uint i=0; i<errors.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanLines(errors[i].second, newLines[i]);
		else if(refStrat_ == LR_SAFE) 
			getFullspanLines(errors[i].second, newLines[i]);
		else if(refStrat_ == LR_ISOTROPIC_FUNC) 
			getStructMeshLines(errors[i].second, newLines[i]);
		// note that this is an excessive loop as it computes the meshlines for ALL elements,
		// but we're only going to use a small part of this.
	}

	/* Do the actual refinement */
	int target_n_functions = ceil(basis_.size()*(1+beta));
	int i=0;
	while( basis_.size() < target_n_functions ) {
		for(uint j=0; j<newLines[i].size(); j++) {
			Meshline *m = newLines[i][j];
			insert_line(!m->is_spanning_u(), m->const_par_, m->start_, m->stop_, refKnotlineMult_);
		}
		i++;
	}

	/* do a posteriori fixes to ensure a proper mesh */
	aPosterioriFixes();

	/* exit cleanly be deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++) 
		for(uint j=0; j<newLines[i].size(); j++) 
			delete newLines[i][j];
}

void LRSplineSurface::aPosterioriFixes()  {
	std::vector<Meshline*> *newLines = NULL;
	int nFunc;
	do {
		nFunc = basis_.size();
		if(doCloseGaps_)
			this->closeGaps(newLines);
		if(maxTjoints_ > 0)
			this->enforceMaxTjoints(newLines);
		if(doAspectRatioFix_)
			this->enforceMaxAspectRatio(newLines);
	} while(nFunc != basis_.size());
}


void LRSplineSurface::closeGaps(std::vector<Meshline*>* newLines) {
	std::vector<double>  start_v;
	std::vector<double>  stop_v ;
	std::vector<double>  const_u  ;

	std::vector<double>  start_u;
	std::vector<double>  stop_u ;
	std::vector<double>  const_v  ;
	for(uint i=0; i<element_.size(); i++) {
		double umin = element_[i]->umin();
		double umax = element_[i]->umax();
		double vmin = element_[i]->vmin();
		double vmax = element_[i]->vmax();
		std::vector<double> left, right, top, bottom;
		for(uint j=0; j<meshline_.size(); j++) {
			Meshline *m = meshline_[j];
			if(      m->span_u_line_ &&
			         m->const_par_ > vmin &&
			         m->const_par_ < vmax) {
				if(m->start_ == umax)
					right.push_back(m->const_par_);
				else if(m->stop_ == umin)
					left.push_back(m->const_par_);
			} else if(!m->span_u_line_ &&
			          m->const_par_ > umin &&
			          m->const_par_ < umax) {
				if(m->start_ == vmax)
					top.push_back(m->const_par_);
				else if(m->stop_ == vmin)
					bottom.push_back(m->const_par_);
			}
		}
		for(uint j=0; j<left.size(); j++)
			for(uint k=0; k<right.size(); k++)
				if(left[j] == right[k]) {
					const_v.push_back(left[j]);
					start_u.push_back(umin);
					stop_u.push_back(umax);
				}
		for(uint j=0; j<top.size(); j++)
			for(uint k=0; k<bottom.size(); k++)
				if(top[j] == bottom[k]) {
					const_u.push_back(top[j]);
					start_v.push_back(vmin);
					stop_v.push_back(vmax);
				}
	}
	Meshline* m;
	for(uint i=0; i<const_u.size(); i++) {
		m = insert_const_u_edge(const_u[i], start_v[i], stop_v[i], refKnotlineMult_);
		if(newLines != NULL)
			newLines->push_back(m->copy());
	}
	for(uint i=0; i<const_v.size(); i++) {
		m = insert_const_v_edge(const_v[i], start_u[i], stop_u[i], refKnotlineMult_);
		if(newLines != NULL)
			newLines->push_back(m->copy());
	}
}

void LRSplineSurface::enforceMaxTjoints(std::vector<Meshline*> *newLines) {
	bool someFix = true;
	while(someFix) {
		someFix = false;
		for(uint i=0; i<element_.size(); i++) {
			double umin = element_[i]->umin();
			double umax = element_[i]->umax();
			double vmin = element_[i]->vmin();
			double vmax = element_[i]->vmax();
			std::vector<double> left, right, top, bottom;
			for(uint j=0; j<meshline_.size(); j++) {
				Meshline *m = meshline_[j];
				if(      m->span_u_line_ &&
						 m->const_par_ > vmin &&
						 m->const_par_ < vmax) {
					if(m->start_ == umax)
						right.push_back(m->const_par_);
					else if(m->stop_ == umin)
						left.push_back(m->const_par_);
				} else if(!m->span_u_line_ &&
						  m->const_par_ > umin &&
						  m->const_par_ < umax) {
					if(m->start_ == vmax)
						top.push_back(m->const_par_);
					else if(m->stop_ == vmin)
						bottom.push_back(m->const_par_);
				}
			}
			Meshline *m;
			double best = DBL_MAX;
			int bi      = -1;
			if(left.size() > (uint) maxTjoints_) {
				for(uint j=0; j<left.size(); j++) {
					if(MY_STUPID_FABS(left[j] - (vmin+vmax)/2) < best) {
						best = MY_STUPID_FABS(left[j] - (vmin+vmax)/2);
						bi = j;
					}
				}
				m = insert_const_v_edge(left[bi], umin, umax, refKnotlineMult_);
				if(newLines != NULL)
					newLines->push_back(m->copy());
				someFix = true;
				continue;
			}
			if(right.size() > (uint) maxTjoints_) {
				for(uint j=0; j<right.size(); j++) {
					if(MY_STUPID_FABS(right[j] - (vmin+vmax)/2) < best) {
						best = MY_STUPID_FABS(right[j] - (vmin+vmax)/2);
						bi = j;
					}
				}
				m = insert_const_v_edge(right[bi], umin, umax, refKnotlineMult_);
				if(newLines != NULL)
					newLines->push_back(m->copy());
				someFix = true;
				continue;
			}
			if(top.size() > (uint) maxTjoints_) {
				for(uint j=0; j<top.size(); j++) {
					if(MY_STUPID_FABS(top[j] - (umin+umax)/2) < best) {
						best = MY_STUPID_FABS(top[j] - (umin+umax)/2);
						bi = j;
					}
				}
				m = insert_const_u_edge(top[bi], vmin, vmax, refKnotlineMult_);
				if(newLines != NULL)
					newLines->push_back(m->copy());
				someFix = true;
				continue;
			}
			if(bottom.size() > (uint) maxTjoints_) {
				for(uint j=0; j<bottom.size(); j++) {
					if(MY_STUPID_FABS(bottom[j] - (umin+umax)/2) < best) {
						best = MY_STUPID_FABS(bottom[j] - (umin+umax)/2);
						bi = j;
					}
				}
				m = insert_const_u_edge(bottom[bi], vmin, vmax, refKnotlineMult_);
				if(newLines != NULL)
					newLines->push_back(m->copy());
				someFix = true;
				continue;
			}
		}
	}
}

void LRSplineSurface::enforceMaxAspectRatio(std::vector<Meshline*>* newLines) {
	bool somethingFixed = true;
	while(somethingFixed) {
		somethingFixed = false;
		for(uint i=0; i<element_.size(); i++) {
			double umin = element_[i]->umin();
			double umax = element_[i]->umax();
			double vmin = element_[i]->vmin();
			double vmax = element_[i]->vmax();
			bool insert_const_u =  umax-umin > maxAspectRatio_*(vmax-vmin);
			bool insert_const_v =  vmax-vmin > maxAspectRatio_*(umax-umin);
			if( insert_const_u || insert_const_v ) {
				std::vector<Meshline*> splitLines; // should always contain exactly one meshline on function return
				if(refStrat_ == LR_MINSPAN) 
					getMinspanLines(i, splitLines);
				else
					getFullspanLines(i, splitLines);

				
				Meshline *m, *msplit;
				msplit = splitLines.front();

				m = insert_line(!msplit->is_spanning_u(), msplit->const_par_, msplit->start_, msplit->stop_, refKnotlineMult_);
				if(newLines != NULL)
					newLines->push_back(m->copy());

				delete msplit;
				somethingFixed = true;
				break;
			}
		}
	}
}

Meshline* LRSplineSurface::insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity) {
	return insert_line(true, u, start_v, stop_v, multiplicity);
}

Meshline* LRSplineSurface::insert_line(bool const_u, double const_par, double start, double stop, int multiplicity) {
	Meshline *newline;
#ifdef TIME_LRSPLINE
	PROFILE("insert_line()");
#endif
	{ // check if the line is an extension or a merging of existing lines
#ifdef TIME_LRSPLINE
	PROFILE("line verification");
#endif
	newline = new Meshline(!const_u, const_par, start, stop, multiplicity);
	newline->type_ = NEWLINE;
	for(uint i=0; i<meshline_.size(); i++) {
		// if newline overlaps any existing ones (may be multiple existing ones)
		// let newline be the entire length of all merged and delete the unused ones

		if(meshline_[i]->is_spanning_u() != const_u && MY_STUPID_FABS(meshline_[i]->const_par_-const_par)<DOUBLE_TOL && 
		   meshline_[i]->start_ <= stop && meshline_[i]->stop_ >= start)  { // meshline_[i] overlaps with newline

			if(meshline_[i]->start_ <= start && 
			   meshline_[i]->stop_  >= stop ) { // newline completely contained in meshline_[i]
			   
				if(meshline_[i]->multiplicity_ < newline->multiplicity_) { // increasing multiplicity
					if(meshline_[i]->start_ == start && 
					   meshline_[i]->stop_  == stop ) { // increasing the mult of the entire line

						// keeping newline, getting rid of the old line
						delete meshline_[i];
						meshline_.erase(meshline_.begin() + i);
						i--;

					} else { // increasing multiplicity of partial line 
						// do nothing. Keep the entire length meshline_[i], and add newline

					}

				} else { // line exist already, do nothing
					delete newline;
					return meshline_[i];
				}
			} else { // newline overlaps meshline_[i]. Keep (and update) newline, delete meshline_[i]

				// update refinement type (for later analysis of linear independence)
				if(newline->type_ == ELONGATION)   // overlaps two existing lines => MERGING
					newline->type_ = MERGING;
				else if(newline->type_ != MERGING) // overlaps one existing line => ELONGATION
					newline->type_ = ELONGATION;

				// update the length of the line with the lowest multiplicity
				if(meshline_[i]->multiplicity_ < newline->multiplicity_) {
					if(meshline_[i]->start_ > start) meshline_[i]->start_ = newline->start_;
					if(meshline_[i]->stop_  < stop ) meshline_[i]->stop_  = newline->stop_;

				} else if(meshline_[i]->multiplicity_ > newline->multiplicity_) {
					if(meshline_[i]->start_ < start) newline->start_ = meshline_[i]->start_;
					if(meshline_[i]->stop_  > stop ) newline->stop_  = meshline_[i]->stop_;

				} else { // for equal mult, we only keep newline and remove the previous line
					if(meshline_[i]->start_ < start) newline->start_ = meshline_[i]->start_;
					if(meshline_[i]->stop_  > stop ) newline->stop_  = meshline_[i]->stop_;

					// keeping newline, getting rid of the old line
					delete meshline_[i];
					meshline_.erase(meshline_.begin() + i);
					i--;
				} 

			}
		}
	}
	}

	HashSet<Basisfunction*> newFuncStp1, newFuncStp2;
	HashSet<Basisfunction*> removeFunc;

	{ // STEP 1: test EVERY function against the NEW meshline
#ifdef TIME_LRSPLINE
	PROFILE("STEP 1");
#endif
	{
#ifdef TIME_LRSPLINE
	PROFILE("S1-basissplit");
#endif
	for(Basisfunction* b : basis_) {
		if(newline->splits(b)) {
			int nKnots;
			if( ((nKnots=newline->nKnotsIn(b)) != newline->multiplicity_) ) {
				removeFunc.insert(b);
				split( const_u, b, const_par, newline->multiplicity_-nKnots, newFuncStp1 );
			}
		}
	}
	for(Basisfunction* b : removeFunc) {
		basis_.erase(b);
		delete b;
	}
	} // end profiler
	{
#ifdef TIME_LRSPLINE
	PROFILE("S1-elementsplit");
#endif
	for(uint i=0; i<element_.size(); i++) {
		if(newline->splits(element_[i]))
			element_.push_back(element_[i]->split(newline->is_spanning_u(), newline->const_par_));
	}
	} // end profiler (elementsplit)
	} // end profiler (step 1)

	{ // STEP 2: test every NEW function against ALL old meshlines
#ifdef TIME_LRSPLINE
	PROFILE("STEP 2");
#endif
	meshline_.push_back(newline);
	while(newFuncStp1.size() > 0) {
		Basisfunction *b = newFuncStp1.pop();
		bool splitMore = false;
		for(Meshline *m : meshline_) {
			if(m->splits(b)) {
				int nKnots = m->nKnotsIn(b);
				if( nKnots != m->multiplicity_ ) {
					splitMore = true;
					split( !m->is_spanning_u(), b, m->const_par_, m->multiplicity_-nKnots, newFuncStp1);
					delete b;
					break;
				}
			}
		}
		if(!splitMore)
			basis_.insert(b);
	}
	} // end profiler (step 2)

	return newline;
}

Meshline* LRSplineSurface::insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity) {
	return insert_line(false, v, start_u, stop_u, multiplicity);
}

void LRSplineSurface::split(bool insert_in_u, Basisfunction* b, double new_knot, int multiplicity, HashSet<Basisfunction*> &newFunctions) {
#ifdef TIME_LRSPLINE
	PROFILE("split()");
#endif

	// create the new functions b1 and b2
	Basisfunction *b1, *b2;
	std::vector<double> knot = (insert_in_u) ? (*b)[0]  : (*b)[1];
	int     p                = (insert_in_u) ? order_u_ : order_v_;
	int     insert_index = 0;
	if(new_knot < knot[0] || knot[p] < new_knot)
		return ;
	while(knot[insert_index] < new_knot)
		insert_index++;
	double alpha1 = (insert_index == p)  ? 1.0 : (new_knot-knot[0])/(knot[p-1]-knot[0]);
	double alpha2 = (insert_index == 1 ) ? 1.0 : (knot[p]-new_knot)/(knot[p]-knot[1]);
	double newKnot[p+2];
	std::copy(knot.begin(), knot.end(), newKnot+1);
	newKnot[0] = new_knot;
	std::sort(newKnot, newKnot + p+2);
	if(insert_in_u) {
		b1 = new Basisfunction(newKnot  , (*b)[1].begin(), b->cp(), b->dim(), order_u_, order_v_, b->w()*alpha1);
		b2 = new Basisfunction(newKnot+1, (*b)[1].begin(), b->cp(), b->dim(), order_u_, order_v_, b->w()*alpha2);
	} else { // insert in v
		b1 = new Basisfunction((*b)[0].begin(), newKnot  , b->cp(), b->dim(), order_u_, order_v_, b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), newKnot+1, b->cp(), b->dim(), order_u_, order_v_, b->w()*alpha2);
	}

	// add any brand new functions and detect their support elements
	HashSet_iterator<Basisfunction*> it = basis_.find(b1);
	if(it != basis_.end()) {
		**it += *b1;
		delete b1;
	} else {
		it = newFunctions.find(b1);
		if(it != newFunctions.end()) {
			**it += *b1;
			delete b1;
		} else {
			updateSupport(b1, b->supportedElementBegin(), b->supportedElementEnd());
			bool recursive_split = (multiplicity > 1) && ( ( insert_in_u && (*b1)[0][order_u_]!=new_knot) ||
			                                               (!insert_in_u && (*b1)[1][order_v_]!=new_knot)  );
			if(recursive_split) {
				split( insert_in_u, b1, new_knot, multiplicity-1, newFunctions);
				delete b1;
			} else {
				newFunctions.insert(b1);
			}
		}
	}
	it = basis_.find(b2);
	if(it != basis_.end()) {
		**it += *b2;
		delete b2;
	} else {
		it = newFunctions.find(b2);
		if(it != newFunctions.end()) {
			**it += *b2;
			delete b2;
		} else {
			updateSupport(b2, b->supportedElementBegin(), b->supportedElementEnd());
			bool recursive_split = (multiplicity > 1) && ( ( insert_in_u && (*b2)[0][0]!=new_knot) ||
			                                               (!insert_in_u && (*b2)[1][0]!=new_knot)  );
			if(recursive_split) {
				split( insert_in_u, b2, new_knot, multiplicity-1, newFunctions);
				delete b2;
			} else {
				newFunctions.insert(b2);
			}
		}
	}

}

void LRSplineSurface::getEdgeFunctions(std::vector<Basisfunction*> &edgeFunctions, parameterEdge edge, int depth) const {
	edgeFunctions.clear();
	for(Basisfunction *b : basis_) {
		switch(edge) {
		case WEST       :
			if((*b)[0][order_u_-depth] == start_u_)
				edgeFunctions.push_back(b);
			break;
		case EAST       :
			if((*b)[0][depth] == end_u_)
				edgeFunctions.push_back(b);
			break;
		case SOUTH      :
			if((*b)[1][order_v_-depth] == start_v_)
				edgeFunctions.push_back(b);
			break;
		case NORTH      :
			if((*b)[1][depth] == end_v_)
				edgeFunctions.push_back(b);
			break;
		case SOUTH_WEST :
			if((*b)[0][order_u_-depth] == start_u_ &&
			   (*b)[1][order_v_-depth] == start_v_)
				edgeFunctions.push_back(b);
			break;
		case SOUTH_EAST :
			if((*b)[0][depth]          == end_u_ &&
			   (*b)[1][order_v_-depth] == start_v_)
				edgeFunctions.push_back(b);
			break;
		case NORTH_WEST :
			if((*b)[0][order_u_-depth] == start_u_ &&
			   (*b)[1][depth]          == end_v_)
				edgeFunctions.push_back(b);
			break;
		case NORTH_EAST :
			if((*b)[0][depth] == end_u_ &&
			   (*b)[1][depth] == end_v_)
				edgeFunctions.push_back(b);
			break;
		default:
			break;
		}
	}
}

void LRSplineSurface::rebuildDimension(int dimvalue) {
	for(Basisfunction *b : basis_)
		b->setDimension(dimvalue);
	dim_ = dimvalue;
}

void LRSplineSurface::getGlobalKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
	getGlobalUniqueKnotVector(knot_u, knot_v);

	// add in duplicates where apropriate
	for(uint i=0; i<knot_u.size(); i++) {
		for(uint j=0; j<meshline_.size(); j++) {
			if(!meshline_[j]->is_spanning_u() && meshline_[j]->const_par_==knot_u[i]) {
				for(int m=1; m<meshline_[j]->multiplicity_; m++) {
					knot_u.insert(knot_u.begin()+i, knot_u[i]);
					i++;
				}
				break;
			}
		}
	}
	for(uint i=0; i<knot_v.size(); i++) {
		for(uint j=0; j<meshline_.size(); j++) {
			if(meshline_[j]->is_spanning_u() && meshline_[j]->const_par_==knot_v[i]) {
				for(int m=1; m<meshline_[j]->multiplicity_; m++) {
					knot_v.insert(knot_v.begin()+i, knot_v[i]);
					i++;
				}
				break;
			}
		}
	}
}

void LRSplineSurface::getGlobalUniqueKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
	knot_u.clear();
	knot_v.clear();
	// create a huge list of all line instances
	for(uint i=0; i<meshline_.size(); i++) {
		if(meshline_[i]->is_spanning_u())
			knot_v.push_back(meshline_[i]->const_par_);
		else
			knot_u.push_back(meshline_[i]->const_par_);
	}

	// sort them and remove duplicates
	std::vector<double>::iterator it;
	std::sort(knot_u.begin(), knot_u.end());
	std::sort(knot_v.begin(), knot_v.end());
	it = std::unique(knot_u.begin(), knot_u.end());
	knot_u.resize(it-knot_u.begin());
	it = std::unique(knot_v.begin(), knot_v.end());
	knot_v.resize(it-knot_v.begin());
}

void LRSplineSurface::updateSupport(Basisfunction *f,
	                                std::vector<Element*>::iterator start,
	                                std::vector<Element*>::iterator end ) {
#ifdef TIME_LRSPLINE
	PROFILE("update support");
#endif
	std::vector<Element*>::iterator it;
	for(it=start; it!=end; it++)
		if(f->addSupport(*it)) // this tests for overlapping as well as updating  
			(*it)->addSupportFunction(f);
}

void LRSplineSurface::updateSupport(Basisfunction *f) {
	updateSupport(f, element_.begin(), element_.end());
}

void LRSplineSurface::generateIDs() const {
	uint i=0;
	for(Basisfunction *b : basis_)
		b->setId(i++);
	for(i=0; i<element_.size(); i++) 
		element_[i]->setId(i);
}

bool LRSplineSurface::isLinearIndepByOverloading(bool verbose) {
	std::vector<Basisfunction*>           overloaded;
	std::vector<int>                      singleElms;
	std::vector<int>                      multipleElms;
	std::vector<int>                      singleBasis;
	std::vector<int>                      multipleBasis;
	// reset overload counts and initialize overloaded basisfunction set
	for(uint i=0; i<element_.size(); i++)
		element_[i]->resetOverloadCount();
	for(Basisfunction *b : basis_)
		if(b->isOverloaded())
			overloaded.push_back(b);
	
	int lastOverloadCount = overloaded.size();
	int iterationCount = 0 ;
	do {
		lastOverloadCount = overloaded.size();
		// reset overload counts and plotting vectors
		for(uint i=0; i<element_.size(); i++)
			element_[i]->resetOverloadCount();
		singleBasis.clear();
		multipleBasis.clear();
		singleElms.clear();
		multipleElms.clear();
		for(uint i=0; i<overloaded.size(); i++)
			for(Element* e : overloaded[i]->support()) 
				e->incrementOverloadCount();
		for(uint i=0; i<element_.size(); i++) {
			if(element_[i]->getOverloadCount() > 1)
				multipleElms.push_back(i);
			if(element_[i]->getOverloadCount() == 1)
				singleElms.push_back(i);
		}
		for(uint i=0; i<overloaded.size(); i++) {
			if(overloaded[i]->getOverloadCount() > 1)
				multipleBasis.push_back(overloaded[i]->getId());
			if(overloaded[i]->getOverloadCount() == 1) {
				singleBasis.push_back(overloaded[i]->getId());
				overloaded.erase(overloaded.begin() + i);
				i--;
			}
		}
		int nOverloadedElms  = singleElms.size() + multipleElms.size();
		if(verbose) {
			std::cout << "Iteration #"<< iterationCount << "\n";
			std::cout << "Overloaded elements  : " << nOverloadedElms
		          	<< " ( " << singleElms.size() << " + "
				  	<< multipleElms.size() << ")\n";
			std::cout << "Overloaded B-splines : " << lastOverloadCount
		          	<< " ( " << overloaded.size() << " + "
				  	<< (lastOverloadCount-overloaded.size()) << ")\n";
			std::cout << "-----------------------------------------------\n\n";

			char filename[256];
			std::ofstream out;
			sprintf(filename, "elements_%03d.eps", iterationCount);
			out.open(filename);
			setElementColor(0.6, 0.6, 0.6);
			writePostscriptElements(out, 2,2,false, &singleElms);
			setElementColor(0.9, 0.3, 0.15);
			writePostscriptElements(out, 2,2,true,  &multipleElms);
			out.close();
		
			sprintf(filename, "functions_%03d.eps", iterationCount);
			out.open(filename);
			setBasisColor(0.6, 0.6, 0.6);
			setSelectedBasisColor(0.9, 0.83, 0.05);
			writePostscriptFunctionSpace(out, &singleBasis, true, false);
			setSelectedBasisColor(1.0, 0.10, 0.00);
			writePostscriptFunctionSpace(out,  &multipleBasis, false, true);
			out.close();
		}


		iterationCount++;
	} while((uint) lastOverloadCount != overloaded.size());

	return overloaded.size() == 0;
}

bool LRSplineSurface::isLinearIndepByMappingMatrix(bool verbose) const {
#ifdef TIME_LRSPLINE
	PROFILE("Linear independent)");
#endif
	// try and figure out this thing by the projection matrix C

	std::vector<double> knots_u, knots_v;
	getGlobalKnotVector(knots_u, knots_v);
	int nmb_bas = basis_.size();
	int n1 = knots_u.size() - order_u_;
	int n2 = knots_v.size() - order_v_;
	int fullDim = n1*n2;
	bool fullVerbose   = fullDim < 30  && nmb_bas < 50;
	bool sparseVerbose = fullDim < 250 && nmb_bas < 100;

	std::vector<std::vector<boost::rational<long long> > > C;  // rational projection matrix 

	// scaling factor to ensure that all knots are integers (assuming all multiplum of smallest knot span)
	double smallKnotU = DBL_MAX;
	double smallKnotV = DBL_MAX;
	for(uint i=0; i<knots_u.size()-1; i++)
		if(knots_u[i+1]-knots_u[i] < smallKnotU && knots_u[i+1] != knots_u[i])
			smallKnotU = knots_u[i+1]-knots_u[i];
	for(uint i=0; i<knots_v.size()-1; i++)
		if(knots_v[i+1]-knots_v[i] < smallKnotV && knots_v[i+1] != knots_v[i])
			smallKnotV = knots_v[i+1]-knots_v[i];

	// HashSet_iterator<Basisfunction*>       it_begin  = basis_.begin();
	// HashSet_iterator<Basisfunction*>       it_end    = basis_.end();
	HashSet_const_iterator<Basisfunction*> cit_begin = basis_.begin();
	HashSet_const_iterator<Basisfunction*> cit_end   = basis_.end();
	HashSet_const_iterator<Basisfunction*> it;
	// for(Basisfunction *b : basis_) {
	for(it=cit_begin; it != cit_end; ++it) {
		Basisfunction *b = *it;
		int startU, startV;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == (*b)[0][0]) break;
		for(int j=0; j<order_u_; j++) {
			if(knots_u[startU] == (*b)[0][j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == (*b)[1][0]) break;
		for(int j=0; j<order_v_; j++) {
			if(knots_v[startV] == (*b)[1][j]) startV--;
			else break;
		}
		startV++;

		std::vector<boost::rational<long long> > rowU(1,1), rowV(1,1);
		int curU = startU+1;
		for(uint j=0; j<locKnotU.size()-1; j++, curU++) {
			if(locKnotU[j+1] != knots_u[curU]) {
				std::vector<boost::rational<long long> > newRowU(rowU.size()+1, boost::rational<long long>(0));
				for(uint k=0; k<rowU.size(); k++) {
					#define U(x) ((long long) (locKnotU[x+k]/smallKnotU + 0.5))
					long long z = (long long) (knots_u[curU] / smallKnotU + 0.5);
					int p = order_u_-1;
					if(z < U(0) || z > U(p+1)) {
						newRowU[k] = rowU[k];
						continue;
					}
					boost::rational<long long> alpha1 = (U(p) <=  z  ) ? 1 : boost::rational<long long>(   z    - U(0),  U(p)  - U(0));
					boost::rational<long long> alpha2 = (z    <= U(1)) ? 1 : boost::rational<long long>( U(p+1) - z   , U(p+1) - U(1));
					newRowU[k]   += rowU[k]*alpha1;
					newRowU[k+1] += rowU[k]*alpha2;
					#undef U
				}
				locKnotU.insert(locKnotU.begin()+(j+1), knots_u[curU]);
				rowU = newRowU;
			}
		}
		int curV = startV+1;
		for(uint j=0; j<locKnotV.size()-1; j++, curV++) {
			if(locKnotV[j+1] != knots_v[curV]) {
				std::vector<boost::rational<long long> > newRowV(rowV.size()+1, boost::rational<long long>(0));
				for(uint k=0; k<rowV.size(); k++) {
					#define V(x) ((long long) (locKnotV[x+k]/smallKnotV + 0.5))
					long long z = (long long) (knots_v[curV] / smallKnotV + 0.5);
					int p = order_v_-1;
					if(z < V(0) || z > V(p+1)) {
						newRowV[k] = rowV[k];
						continue;
					}
					boost::rational<long long> alpha1 = (V(p) <=  z  ) ? 1 : boost::rational<long long>(   z    - V(0),  V(p)  - V(0));
					boost::rational<long long> alpha2 = (z    <= V(1)) ? 1 : boost::rational<long long>( V(p+1) - z   , V(p+1) - V(1));
					newRowV[k]   += rowV[k]*alpha1;
					newRowV[k+1] += rowV[k]*alpha2;
					#undef V
				}
				locKnotV.insert(locKnotV.begin()+(j+1), knots_v[curV]);
				rowV= newRowV;
			}
		}
		std::vector<boost::rational<long long> > totalRow(fullDim, boost::rational<long long>(0));
		for(uint i1=0; i1<rowU.size(); i1++)
			for(uint i2=0; i2<rowV.size(); i2++)
				totalRow[(startV+i2)*n1 + (startU+i1)] = rowV[i2]*rowU[i1];

		C.push_back(totalRow);
	}

	if(verbose && sparseVerbose) {
		for(uint i=0; i<C.size(); i++) {
			std::cout << "|";
			for(uint j=0; j<C[i].size(); j++)
				if(C[i][j]==0) {
					if(fullVerbose)
						std::cout << "\t";
					else if(sparseVerbose)
						std::cout << " ";
				} else {
					if(fullVerbose)
						std::cout << C[i][j] << "\t";
					else if(sparseVerbose)
						std::cout << "x";
				}
			std::cout << "|\n";
		}
		std::cout << std::endl;
	}
	
	// gauss-jordan elimination
	int zeroColumns = 0;
	for(uint i=0; i<C.size() && i+zeroColumns<C[i].size(); i++) {
		boost::rational<long long> maxPivot = 0;
		int maxI = -1;
		for(uint j=i; j<C.size(); j++) {
			if(abs(C[j][i+zeroColumns]) > maxPivot) {
				maxPivot = abs(C[j][i+zeroColumns]);
				maxI = j;
			}
		}
		if(maxI == -1) {
			i--;
			zeroColumns++;
			continue;
		}
		std::vector<boost::rational<long long> > tmp = C[i];
		C[i] = C[maxI];
		C[maxI] = tmp;
		for(uint j=i+1; j<C.size(); j++) {
			// if(j==i) continue;
			boost::rational<long long> scale =  C[j][i+zeroColumns] / C[i][i+zeroColumns];
			if(scale != 0) {
				// std::cout << scale << " x row(" << i << ") added to row " << j << std::endl;
				for(uint k=i+zeroColumns; k<C[j].size(); k++)
					C[j][k] -= C[i][k] * scale;
			}
		}
	}

	if(verbose && sparseVerbose) {
		for(uint i=0; i<C.size(); i++) {
			std::cout << "|";
			for(uint j=0; j<C[i].size(); j++)
				if(C[i][j]==0) {
					if(fullVerbose)
						std::cout << "\t";
					else if(sparseVerbose)
				 	std::cout << " ";
				} else {
					if(fullVerbose)
						std::cout << C[i][j] << "\t";
					else if(sparseVerbose)
						std::cout << "x";
				}
			std::cout << "|\n";
		}
		std::cout << std::endl;
	}

	int rank = (nmb_bas < n1*n2-zeroColumns) ? nmb_bas : n1*n2-zeroColumns;
	if(verbose) {
		std::cout << "Matrix size : " << nmb_bas << " x " << n1*n2 << std::endl;
		std::cout << "Matrix rank : " << rank << std::endl;
	}
	
	return rank == nmb_bas;
}

void LRSplineSurface::getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const {
#ifdef TIME_LRSPLINE
	PROFILE("Linear independent)");
#endif
	// try and figure out this thing by the projection matrix C

	std::vector<double> knots_u, knots_v;
	getGlobalKnotVector(knots_u, knots_v);
	int nmb_bas = basis_.size();
	int n1 = knots_u.size() - order_u_;
	int n2 = knots_v.size() - order_v_;
	int fullDim = n1*n2;

	std::vector<std::vector<boost::rational<long long> > > Ct;  // rational projection matrix (transpose of this)

	// scaling factor to ensure that all knots are integers (assuming all multiplum of smallest knot span)
	double smallKnotU = DBL_MAX;
	double smallKnotV = DBL_MAX;
	for(uint i=0; i<knots_u.size()-1; i++)
		if(knots_u[i+1]-knots_u[i] < smallKnotU && knots_u[i+1] != knots_u[i])
			smallKnotU = knots_u[i+1]-knots_u[i];
	for(uint i=0; i<knots_v.size()-1; i++)
		if(knots_v[i+1]-knots_v[i] < smallKnotV && knots_v[i+1] != knots_v[i])
			smallKnotV = knots_v[i+1]-knots_v[i];

	// initialize C transpose with zeros
	for(int i=0; i<fullDim; i++)
		Ct.push_back(std::vector<boost::rational<long long> >(nmb_bas, boost::rational<long long>(0)));

	int basCount = 0;
	for(Basisfunction *b : basis_)  {
		int startU, startV;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == (*(b))[0][0]) break;
		for(int j=0; j<order_u_; j++) {
			if(knots_u[startU] == (*(b))[0][j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == (*(b))[1][0]) break;
		for(int j=0; j<order_v_; j++) {
			if(knots_v[startV] == (*(b))[1][j]) startV--;
			else break;
		}
		startV++;

		std::vector<boost::rational<long long> > rowU(1,1), rowV(1,1);
		int curU = startU+1;
		for(uint j=0; j<locKnotU.size()-1; j++, curU++) {
			if(locKnotU[j+1] != knots_u[curU]) {
				std::vector<boost::rational<long long> > newRowU(rowU.size()+1, boost::rational<long long>(0));
				for(uint k=0; k<rowU.size(); k++) {
					#define U(x) ((long long) (locKnotU[x+k]/smallKnotU + 0.5))
					long long z = (long long) (knots_u[curU] / smallKnotU + 0.5);
					int p = order_u_-1;
					if(z < U(0) || z > U(p+1)) {
						newRowU[k] = rowU[k];
						continue;
					}
					boost::rational<long long> alpha1 = (U(p) <=  z  ) ? 1 : boost::rational<long long>(   z    - U(0),  U(p)  - U(0));
					boost::rational<long long> alpha2 = (z    <= U(1)) ? 1 : boost::rational<long long>( U(p+1) - z   , U(p+1) - U(1));
					newRowU[k]   += rowU[k]*alpha1;
					newRowU[k+1] += rowU[k]*alpha2;
					#undef U
				}
				locKnotU.insert(locKnotU.begin()+(j+1), knots_u[curU]);
				rowU = newRowU;
			}
		}
		int curV = startV+1;
		for(uint j=0; j<locKnotV.size()-1; j++, curV++) {
			if(locKnotV[j+1] != knots_v[curV]) {
				std::vector<boost::rational<long long> > newRowV(rowV.size()+1, boost::rational<long long>(0));
				for(uint k=0; k<rowV.size(); k++) {
					#define V(x) ((long long) (locKnotV[x+k]/smallKnotV + 0.5))
					long long z = (long long) (knots_v[curV] / smallKnotV + 0.5);
					int p = order_v_-1;
					if(z < V(0) || z > V(p+1)) {
						newRowV[k] = rowV[k];
						continue;
					}
					boost::rational<long long> alpha1 = (V(p) <=  z  ) ? 1 : boost::rational<long long>(   z    - V(0),  V(p)  - V(0));
					boost::rational<long long> alpha2 = (z    <= V(1)) ? 1 : boost::rational<long long>( V(p+1) - z   , V(p+1) - V(1));
					newRowV[k]   += rowV[k]*alpha1;
					newRowV[k+1] += rowV[k]*alpha2;
					#undef V
				}
				locKnotV.insert(locKnotV.begin()+(j+1), knots_v[curV]);
				rowV= newRowV;
			}
		}
		for(uint i1=0; i1<rowU.size(); i1++)
			for(uint i2=0; i2<rowV.size(); i2++)
				Ct[(startV+i2)*n1 + (startU+i1)][basCount] = rowV[i2]*rowU[i1];
		basCount++;
	}


	std::ofstream out;
	out.open("Ct.m");
	out << "A = [";
	for(uint i=0; i<Ct.size(); i++) {
		for(uint j=0; j<Ct[i].size(); j++)
			if(Ct[i][j] > 0)
				out << i << " " << j << " " << Ct[i][j] << ";\n";
	}
	out << "];" << std::endl;
	out.close();
	// gauss-jordan elimination
	int zeroColumns = 0;
	for(uint i=0; i<Ct.size() && i+zeroColumns<Ct[i].size(); i++) {
		boost::rational<long long> maxPivot(0);
		int maxI = -1;
		for(uint j=i; j<Ct.size(); j++) {
			if(abs(Ct[j][i+zeroColumns]) > maxPivot) {
				maxPivot = abs(Ct[j][i+zeroColumns]);
				maxI = j;
			}
		}
		if(maxI == -1) {
			i--;
			zeroColumns++;
			continue;
		}
		std::vector<boost::rational<long long> > tmp = Ct[i];
		Ct[i]    = Ct[maxI];
		Ct[maxI] = tmp;
		boost::rational<long long> leading = Ct[i][i+zeroColumns];
		for(uint j=i+zeroColumns; j<Ct[i].size(); j++)
			Ct[i][j] /= leading;

		for(uint j=0; j<Ct.size(); j++) {
			if(i==j) continue;
			boost::rational<long long> scale =  Ct[j][i+zeroColumns] / Ct[i][i+zeroColumns];
			if(scale != 0) {
				for(uint k=i+zeroColumns; k<Ct[j].size(); k++)
					Ct[j][k] -= Ct[i][k] * scale;
			}
		}
	}

	// figure out the null space
	nullspace.clear();
	for(int i=0; i<zeroColumns; i++) {
		std::vector<boost::rational<long long> > newRow(nmb_bas, boost::rational<long long>(0));
		newRow[nmb_bas-(i+1)] = 1;
		nullspace.push_back(newRow);
	}
	for(int i=0; i<zeroColumns; i++)
		for(int j=0; j<nmb_bas-zeroColumns; j++)
			nullspace[i][j] = -Ct[j][nmb_bas-(i+1)]; // no need to divide by C[j][j] since it's 1

}

bool LRSplineSurface::isLinearIndepByFloatingPointMappingMatrix(bool verbose) const {
#ifdef TIME_LRSPLINE
	PROFILE("Linear independent)");
#endif
	// try and figure out this thing by the projection matrix C

	std::vector<double> knots_u, knots_v;
	getGlobalKnotVector(knots_u, knots_v);
	int nmb_bas = basis_.size();
	int n1 = knots_u.size() - order_u_;
	int n2 = knots_v.size() - order_v_;
	int fullDim = n1*n2;
	bool fullVerbose   = fullDim < 30  && nmb_bas < 50;
	bool sparseVerbose = fullDim < 250 && nmb_bas < 100;

	std::vector<std::vector<double > > C;  // rational projection matrix 

	for(Basisfunction *b : basis_)  {
		int startU, startV;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == (*b)[0][0]) break;
		for(int j=0; j<order_u_; j++) {
			if(knots_u[startU] == (*b)[0][j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == (*b)[1][0]) break;
		for(int j=0; j<order_v_; j++) {
			if(knots_v[startV] == (*b)[1][j]) startV--;
			else break;
		}
		startV++;

		std::vector<double > rowU(1,1), rowV(1,1);
		int curU = startU+1;
		for(uint j=0; j<locKnotU.size()-1; j++, curU++) {
			if(locKnotU[j+1] != knots_u[curU]) {
				std::vector<double > newRowU(rowU.size()+1, 0);
				for(uint k=0; k<rowU.size(); k++) {
					#define U(x) (locKnotU[x+k])
					double z = knots_u[curU] ;
					int p = order_u_-1;
					if(z < U(0) || z > U(p+1)) {
						newRowU[k] = rowU[k];
						continue;
					}
					double alpha1 = (U(p) <=  z  ) ? 1 : (   z    - U(0))/(  U(p)  - U(0));
					double alpha2 = (z    <= U(1)) ? 1 : ( U(p+1) - z   )/( U(p+1) - U(1));
					newRowU[k]   += rowU[k]*alpha1;
					newRowU[k+1] += rowU[k]*alpha2;
					#undef U
				}
				locKnotU.insert(locKnotU.begin()+(j+1), knots_u[curU]);
				rowU = newRowU;
			}
		}
		int curV = startV+1;
		for(uint j=0; j<locKnotV.size()-1; j++, curV++) {
			if(locKnotV[j+1] != knots_v[curV]) {
				std::vector<double > newRowV(rowV.size()+1, 0);
				for(uint k=0; k<rowV.size(); k++) {
					#define V(x) (locKnotV[x+k])
					double z = knots_v[curV];
					int p = order_v_-1;
					if(z < V(0) || z > V(p+1)) {
						newRowV[k] = rowV[k];
						continue;
					}
					double alpha1 = (V(p) <=  z  ) ? 1 : (   z    - V(0))/(  V(p)  - V(0));
					double alpha2 = (z    <= V(1)) ? 1 : ( V(p+1) - z   )/( V(p+1) - V(1));
					newRowV[k]   += rowV[k]*alpha1;
					newRowV[k+1] += rowV[k]*alpha2;
					#undef V
				}
				locKnotV.insert(locKnotV.begin()+(j+1), knots_v[curV]);
				rowV= newRowV;
			}
		}
		std::vector<double > totalRow(fullDim, 0.0);
		for(uint i1=0; i1<rowU.size(); i1++)
			for(uint i2=0; i2<rowV.size(); i2++)
				totalRow[(startV+i2)*n1 + (startU+i1)] = rowV[i2]*rowU[i1];

		C.push_back(totalRow);
	}

	if(verbose && sparseVerbose) {
		for(uint i=0; i<C.size(); i++) {
			std::cout << "|";
			for(uint j=0; j<C[i].size(); j++)
				if(C[i][j]==0) {
					if(fullVerbose)
						std::cout << "\t";
					else if(sparseVerbose)
						std::cout << " ";
				} else {
					if(fullVerbose)
						std::cout << C[i][j] << "\t";
					else if(sparseVerbose)
						std::cout << "x";
				}
			std::cout << "|\n";
		}
		std::cout << std::endl;
	}

	// gauss-jordan elimination
	int zeroColumns = 0;
	for(uint i=0; i<C.size() && i+zeroColumns<C[i].size(); i++) {
		double maxPivot = 0;
		int maxI = -1;
		for(uint j=i; j<C.size(); j++) {
			if(MY_STUPID_FABS(C[j][i+zeroColumns]) > maxPivot) {
				maxPivot = MY_STUPID_FABS(C[j][i+zeroColumns]);
				maxI = j;
			}
		}
		if(maxPivot > 0 && maxPivot < 1e-10) {
			C[maxI][i+zeroColumns] = 0.0;
			maxI = -1;
		}
		if(maxI == -1) {
			i--;
			zeroColumns++;
			continue;
		}
		std::vector<double> tmp = C[i];
		C[i] = C[maxI];
		C[maxI] = tmp;
		for(uint j=i+1; j<C.size(); j++) {
			double scale =  C[j][i+zeroColumns] / C[i][i+zeroColumns];
			if(scale != 0) {
				for(uint k=i+zeroColumns; k<C[j].size(); k++)
					C[j][k] -= C[i][k] * scale;
			}
		}
	}

	if(verbose && sparseVerbose) {
		for(uint i=0; i<C.size(); i++) {
			std::cout << "|";
			for(uint j=0; j<C[i].size(); j++)
				if(C[i][j]==0) {
					if(fullVerbose)
						std::cout << "\t";
					else if(sparseVerbose)
				 	std::cout << " ";
				} else {
					if(fullVerbose)
						std::cout << C[i][j] << "\t";
					else if(sparseVerbose)
						std::cout << "x";
				}
			std::cout << "|\n";
		}
		std::cout << std::endl;
	}

	int rank = (nmb_bas < n1*n2-zeroColumns) ? nmb_bas : n1*n2-zeroColumns;
	if(verbose) {
		std::cout << "Matrix size : " << nmb_bas << " x " << n1*n2 << std::endl;
		std::cout << "Matrix rank : " << rank << std::endl;
	}
	
	return rank == nmb_bas;
}

double LRSplineSurface::makeIntegerKnots() {
	// find the smallest knot interval
	double smallKnotU = DBL_MAX;
	double smallKnotV = DBL_MAX;
	std::vector<double> knots_u, knots_v;
	getGlobalKnotVector(knots_u, knots_v);
	for(uint i=0; i<knots_u.size()-1; i++)
		if(knots_u[i+1]-knots_u[i] < smallKnotU && knots_u[i+1] != knots_u[i])
			smallKnotU = knots_u[i+1]-knots_u[i];
	for(uint i=0; i<knots_v.size()-1; i++)
		if(knots_v[i+1]-knots_v[i] < smallKnotV && knots_v[i+1] != knots_v[i])
			smallKnotV = knots_v[i+1]-knots_v[i];
	double scale = (smallKnotU<smallKnotV) ? smallKnotU : smallKnotV;

	// scale all meshline values according to this
	Meshline *m;
	for(uint i=0; i<meshline_.size(); i++) {
		m = meshline_[i];
		m->const_par_ = floor(m->const_par_/scale + 0.5);
		m->start_     = floor(m->start_    /scale + 0.5);
		m->stop_      = floor(m->stop_     /scale + 0.5);
	}

	// scale all element values 
	Element *e;
	for(uint i=0; i<element_.size(); i++) {
		e = element_[i];
		e->setUmin( floor(e->umin() /scale + 0.5));
		e->setVmin( floor(e->vmin() /scale + 0.5));
		e->setUmax( floor(e->umax() /scale + 0.5));
		e->setVmax( floor(e->vmax() /scale + 0.5));
	}

	// scale all basis functions values
	for(Basisfunction *b : basis_) {
		for(int j=0; j<order_u_+1; j++)
			(*b)[0][j] = floor((*b)[0][j]/scale + 0.5);
		for(int j=0; j<order_v_+1; j++)
			(*b)[1][j] = floor((*b)[1][j]/scale + 0.5);
	}
	
	// scale all LRSplineSurface values
	start_u_ = floor(start_u_/scale + 0.5);
	start_v_ = floor(start_v_/scale + 0.5);
	end_u_   = floor(end_u_  /scale + 0.5);
	end_v_   = floor(end_v_  /scale + 0.5);

	return scale;
}

void LRSplineSurface::getSupportElements(std::vector<int> &result, const std::vector<int> &basisfunctions) const  {
	result.clear();
	std::set<int> tmp;
	std::vector<Element*>::iterator it;
	for(Basisfunction *b : basis_) {
		for(it=b->supportedElementBegin(); it != b->supportedElementEnd(); it++)
			tmp.insert((**it).getId());
	}
	result.resize(tmp.size());
	result.insert(result.begin(), tmp.begin(), tmp.end());
}

void LRSplineSurface::getDiagonalElements(std::vector<int> &result) const  {
	result.clear();
	for(uint i=0; i<element_.size(); i++) 
		if(element_[i]->umin() == element_[i]->vmin())
			result.push_back(i);
}

void LRSplineSurface::getDiagonalBasisfunctions(std::vector<int> &result) const  {
	result.clear();
	int i = 0;
	for(Basisfunction *b : basis_) {
		bool isDiag = true;
		for(int j=0; j<order_u_+1; j++)
			if((*b)[0][j] != (*b)[1][j])
				isDiag = false;
		if(isDiag)
			result.push_back(i);
		i++;
	}
}

void LRSplineSurface::setElementColor(double r, double g, double b)  {
	element_red   = r;
	element_green = g;
	element_blue  = b;
}

void LRSplineSurface::setBasisColor(double r, double g, double b)  {
	basis_red   = r;
	basis_green = g;
	basis_blue  = b;
}

void LRSplineSurface::setSelectedBasisColor(double r, double g, double b) { 
	selected_basis_red   = r;
	selected_basis_green = g;
	selected_basis_blue  = b;
}

void LRSplineSurface::read(std::istream &is) {
	start_u_ =  DBL_MAX;
	end_u_   = -DBL_MAX;
	start_v_ =  DBL_MAX;
	end_v_   = -DBL_MAX;

	// first get rid of comments and spaces
	ws(is); 
	char firstChar;
	char buffer[1024];
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	// read actual parameters
	int nBasis, nElements, nMeshlines;
	is >> order_u_;    ws(is);
	is >> order_v_;    ws(is);
	is >> nBasis;      ws(is);
	is >> nMeshlines;  ws(is);
	is >> nElements;   ws(is);
	is >> dim_;        ws(is);
	is >> rational_;   ws(is);
	
	meshline_.resize(nMeshlines);
	element_.resize(nElements);
	std::vector<Basisfunction*> basisVector(nBasis);

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	// read all basisfunctions
	for(int i=0; i<nBasis; i++) {
		Basisfunction *b = new Basisfunction(dim_, order_u_, order_v_);
		b->read(is);
		basis_.insert(b);
		basisVector[i] = b;
	}

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	for(int i=0; i<nMeshlines; i++) {
		meshline_[i] = new Meshline();
		meshline_[i]->read(is);
	}

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	// read elements and calculate patch boundaries
	for(int i=0; i<nElements; i++) {
		element_[i] = new Element();
		element_[i]->read(is);
		element_[i]->updateBasisPointers(basisVector);
		start_u_ = (element_[i]->umin() < start_u_) ? element_[i]->umin() : start_u_;
		end_u_   = (element_[i]->umax() > end_u_  ) ? element_[i]->umax() : end_u_  ;
		start_v_ = (element_[i]->vmin() < start_v_) ? element_[i]->vmin() : start_v_;
		end_v_   = (element_[i]->vmax() > end_v_  ) ? element_[i]->vmax() : end_v_  ;
	}
}

void LRSplineSurface::write(std::ostream &os) const {
	generateIDs();
	os << std::setprecision(16);
	os << "# LRSPLINE SURFACE\n";
	os << "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n\t";
	os << order_u_ << "\t";
	os << order_v_ << "\t";
	os << basis_.size() << "\t";
	os << meshline_.size() << "\t";
	os << element_.size() << "\t";
	os << dim_ << "\t";
	os << rational_ << "\n";

	os << "# Basis functions:\n";
	for(Basisfunction* b : basis_) 
		os << *b << std::endl;
	os << "# Mesh lines:\n";
	for(Meshline* m : meshline_)  
		os << *m << std::endl;
	os << "# Elements:\n";
	for(Element* e : element_)
		os << *e << std::endl;
}

void LRSplineSurface::writePostscriptMesh(std::ostream &out, bool close, std::vector<int> *colorElements) const {
#ifdef TIME_LRSPLINE
	PROFILE("Write EPS");
#endif
	std::vector<double> knot_u, knot_v;
	getGlobalUniqueKnotVector(knot_u, knot_v);
	double min_span_u = knot_u[1] - knot_u[0];
	double min_span_v = knot_v[1] - knot_v[0];
	for(uint i=1; i<knot_u.size()-1; i++)
		min_span_u = (min_span_u<knot_u[i+1]-knot_u[i]) ? min_span_u : knot_u[i+1]-knot_u[i];
	for(uint i=1; i<knot_v.size()-1; i++)
		min_span_v = (min_span_v<knot_v[i+1]-knot_v[i]) ? min_span_v : knot_v[i+1]-knot_v[i];

	// get date
	time_t t = time(0);
	tm* lt = localtime(&t);
	char date[11];
	sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

	// get bounding box
	int dx = end_u_ - start_u_;
	int dy = end_v_ - start_v_;
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;
	// set the duplicate-knot-line (dkl) display width
	double dkl_range = (min_span_u>min_span_v) ? min_span_v*scale/6.0 : min_span_u*scale/6.0; 
	int xmin = (start_u_ - dx/30.0)*scale;
	int ymin = (start_v_ - dy/30.0)*scale;
	int xmax = (end_u_   + dx/30.0)*scale + dkl_range;
	int ymax = (end_v_   + dy/30.0)*scale + dkl_range;

	// print eps header
	out << "%!PS-Adobe-3.0 EPSF-3.0\n";
	out << "%%Creator: LRSplineHelpers.cpp object\n";
	out << "%%Title: LR-spline parameter domain\n";
	out << "%%CreationDate: " << date << std::endl;
	out << "%%Origin: 0 0\n";
	out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;

	// Fill diagonal elements when refining
	if(colorElements != NULL) {
		out << element_red   << " ";
		out << element_green << " ";
		out << element_blue  << " ";
		out << "setrgbcolor \n";
		for(uint i=0; i<colorElements->size(); i++) {
			Element* e = element_[colorElements->at(i)];
			out << "newpath\n";
			out <<  e->umin()*scale << " " << e->vmin()*scale << " moveto\n";
			out <<  e->umax()*scale << " " << e->vmin()*scale << " lineto\n";
			out <<  e->umax()*scale << " " << e->vmax()*scale << " lineto\n";
			out <<  e->umin()*scale << " " << e->vmax()*scale << " lineto\n";
			out << "closepath\n";
			out << "fill\n";
		}
	}

	out << "0 setgray\n";
	out << "1 setlinewidth\n";
	for(uint i=0; i<meshline_.size(); i++) {
		out << "newpath\n";
		double dm = (meshline_[i]->multiplicity_==1) ? 0 : dkl_range/(meshline_[i]->multiplicity_-1);
		for(int m=0; m<meshline_[i]->multiplicity_; m++) {
			if(meshline_[i]->is_spanning_u()) {
				out << meshline_[i]->start_*scale << " " << meshline_[i]->const_par_*scale + dm*m << " moveto\n";
				if(meshline_[i]->stop_ == end_u_)
					out << meshline_[i]->stop_*scale+dkl_range << " " << meshline_[i]->const_par_*scale + dm*m << " lineto\n";
				else
					out << meshline_[i]->stop_*scale << " " << meshline_[i]->const_par_*scale + dm*m << " lineto\n";
			} else {
				out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->start_*scale << " moveto\n";
				if(meshline_[i]->stop_ == end_v_)
					out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->stop_ *scale+dkl_range << " lineto\n";
				else
					out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->stop_ *scale << " lineto\n";
			}
		}
		out << "stroke\n";
	}

	if(close)
		out << "%%EOF\n";
}

void LRSplineSurface::writePostscriptElements(std::ostream &out, int nu, int nv, bool close, std::vector<int> *colorElements) const {
#ifdef TIME_LRSPLINE
	PROFILE("Write EPS");
#endif

	// get date
	time_t t = time(0);
	tm* lt = localtime(&t);
	char date[11];
	sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

	// get bounding box (max/min of the control points)
	double x[2];
	double y[2];
	x[0] = 1e7;
	x[1] = -1e7;
	y[0] = 1e7;
	y[1] = -1e7;
	for(Basisfunction *b : basis_) {
		std::vector<double>::const_iterator cp = b->cp();
		x[0] = (cp[0] < x[0]) ? cp[0] : x[0];
		x[1] = (cp[0] > x[1]) ? cp[0] : x[1];
		y[0] = (cp[1] < y[0]) ? cp[1] : y[0];
		y[1] = (cp[1] > y[1]) ? cp[1] : y[1];
	}

	double dx = x[1]-x[0];
	double dy = y[1]-y[0];
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;

	int xmin = (x[0] - dx/20.0)*scale;
	int ymin = (y[0] - dy/20.0)*scale;
	int xmax = (x[1] + dx/20.0)*scale;
	int ymax = (y[1] + dy/20.0)*scale;

	// print eps header
	out << "%!PS-Adobe-3.0 EPSF-3.0\n";
	out << "%%Creator: LRSplineHelpers.cpp object\n";
	out << "%%Title: LR-spline physical domain\n";
	out << "%%CreationDate: " << date << std::endl;
	out << "%%Origin: 0 0\n";
	out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;
	
	// Fill diagonal elements when refining
	if(colorElements != NULL) {
		out << element_red   << " ";
		out << element_green << " ";
		out << element_blue  << " ";
		out << "setrgbcolor \n";
		for(uint iEl=0; iEl<element_.size(); iEl++) {
			bool doColor = false;
			for(uint j=0; j<colorElements->size(); j++)
				if(iEl == (uint) colorElements->at(j))
					doColor = true;
			if(doColor) {
				double umin = element_[iEl]->umin();
				double umax = element_[iEl]->umax();
				double vmin = element_[iEl]->vmin();
				double vmax = element_[iEl]->vmax();

				Go::Point pt;
				point(pt, umin, vmin, iEl);
				out << "newpath\n";
				out <<  pt[0]*scale << " " << pt[1]*scale << " moveto\n";
				for(int i=1; i<nu; i++) { // SOUTH side
					double u = umin + (umax-umin)*i/(nu-1);
					double v = vmin;
					point(pt, u, v, iEl, false, true);
					out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
				}
				for(int i=1; i<nv; i++) { // EAST side
					double u = umax;
					double v = vmin + (vmax-vmin)*i/(nv-1);
					point(pt, u, v, iEl, false, false);
					out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
				}
				for(int i=nu-1; i-->0; ) { // NORTH side
					double u = umin + (umax-umin)*i/(nu-1);
					double v = vmax;
					point(pt, u, v, iEl, true, false);
					out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
				}
				for(int i=nv-1; i-->1; ) { // WEST side
					double u = umin;
					double v = vmin + (vmax-vmin)*i/(nv-1);
					point(pt, u, v, iEl, true, false);
					out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
				}
				out << "closepath\n";
				out << "fill\n";
			}
		}
	}

	out << "0 setgray\n";
	out << "1 setlinewidth\n";

	for(uint iEl=0; iEl<element_.size(); iEl++) {
		double umin = element_[iEl]->umin();
		double umax = element_[iEl]->umax();
		double vmin = element_[iEl]->vmin();
		double vmax = element_[iEl]->vmax();


/*       DIAGONAL TEST REFINEMENT
 *       used to color the diagonal elements (only 2x2 evaluation points)
 *       maybe dust off this and color arbitrary elements at a later point
 */
/*		point(corner[0], element_[iEl]->umin(), element_[iEl]->vmin(), iEl);
		point(corner[1], element_[iEl]->umin(), element_[iEl]->vmax(), iEl);
		point(corner[2], element_[iEl]->umax(), element_[iEl]->vmax(), iEl);
		point(corner[3], element_[iEl]->umax(), element_[iEl]->vmin(), iEl);
		
		if(colorDiag && element_[iEl]->umin() == element_[iEl]->vmin()) {
			out << ".5 setgray\n";
			if(element_[iEl]->umin() == element_[iEl]->vmin()) {
				out << "newpath\n";
				out <<  corner[0][0]*scale << " " << corner[0][1]*scale << " moveto\n";
				out <<  corner[1][0]*scale << " " << corner[1][1]*scale << " lineto\n";
				out <<  corner[2][0]*scale << " " << corner[2][1]*scale << " lineto\n";
				out <<  corner[3][0]*scale << " " << corner[3][1]*scale << " lineto\n";
				out << "closepath\n";
				out << "fill\n";
			}
			out << "0 setgray\n";
		}
*/

		Go::Point pt;
		point(pt, umin, vmin, iEl);
		out << "newpath\n";
		out <<  pt[0]*scale << " " << pt[1]*scale << " moveto\n";
		for(int i=1; i<nu; i++) { // SOUTH side
			double u = umin + (umax-umin)*i/(nu-1);
			double v = vmin;
			point(pt, u, v, iEl, false, true);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=1; i<nv; i++) { // EAST side
			double u = umax;
			double v = vmin + (vmax-vmin)*i/(nv-1);
			point(pt, u, v, iEl, false, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=nu-1; i-->0; ) { // NORTH side
			double u = umin + (umax-umin)*i/(nu-1);
			double v = vmax;
			point(pt, u, v, iEl, true, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=nv-1; i-->1; ) { // WEST side
			double u = umin;
			double v = vmin + (vmax-vmin)*i/(nv-1);
			point(pt, u, v, iEl, true, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		out << "closepath\n";
		out << "stroke\n";
	}

	if(close)
		out << "%%EOF\n";
}

void LRSplineSurface::writePostscriptMeshWithControlPoints(std::ostream &out, int nu, int nv) const {
	writePostscriptElements(out, nu, nv, false);

	// get bounding box (max/min of the control points)
	double x[2];
	double y[2];
	x[0] = 1e7;
	x[1] = -1e7;
	y[0] = 1e7;
	y[1] = -1e7;
	for(Basisfunction *b : basis_) {
		std::vector<double>::const_iterator cp = b->cp();
		x[0] = (cp[0] < x[0]) ? cp[0] : x[0];
		x[1] = (cp[0] > x[1]) ? cp[0] : x[1];
		y[0] = (cp[1] < y[0]) ? cp[1] : y[0];
		y[1] = (cp[1] > y[1]) ? cp[1] : y[1];
	}

	double dx = x[1]-x[0];
	double dy = y[1]-y[0];
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;

	double circleSize = 15.0;
	
	// create the ellipse function
	out << "/ellipse {\n";
	out << "/endangle exch def\n";
	out << "/startangle exch def\n";
	out << "/yrad exch def\n";
	out << "/xrad exch def\n";
	out << "/y exch def\n";
	out << "/x exch def\n";
	out << "/savematrix matrix currentmatrix def\n";
	out << "x y translate\n";
	out << "xrad yrad scale\n";
	out << "0 0 1 startangle endangle arc\n";
	out << "savematrix setmatrix\n";
	out << "} def\n";

	// load the font to use
	out << "/Times-Roman findfont\n";
	out << "24 scalefont\n";
	out << "setfont\n";

	int i=-1;
	for(Basisfunction *b : basis_) {
		i++;
		double cp_x = b->cp(0);
		double cp_y = b->cp(1);
		// move C^{-1} text on internal functions
		int textX   = ((*b)[0][1] == (*b)[0][order_u_]) ? -2 : 1;
		int textY   = ((*b)[1][1] == (*b)[1][order_v_]) ? -2 : 1;
		// move text on edge functions
		if((*b)[0][1] == end_u_)
			textX = 1;
		else if((*b)[0][order_u_-1] == start_u_)
			textX = -2;
		if((*b)[1][1] == end_v_)
			textY = 1;
		else if((*b)[1][order_v_-1] == start_v_)
			textY = -2;

		out << "newpath\n";
		out << "0.45 0.45 0.45 setrgbcolor \n";
		out << cp_x*scale << " " << cp_y*scale << " " << circleSize << " " << circleSize << " 0 360 ellipse\n";
		out << "closepath fill\n";
		out << "0 setgray\n";
		out << cp_x*scale << " " << cp_y*scale << " " << circleSize << " " << circleSize << " 0 360 ellipse\n";
		out << "closepath stroke\n";
		out << "\n";
		out << "newpath\n";
		out << cp_x*scale + textX*circleSize << " " << cp_y*scale + textY*circleSize << " moveto\n";
		out << "(" << i << ") show\n";
		out << "\n";
	}
	out << "%%EOF\n";
}

void LRSplineSurface::writePostscriptFunctionSpace(std::ostream &out, std::vector<int> *colorBasis, bool drawAll, bool close) const {
	if(drawAll)
		writePostscriptMesh(out, false);

	int dx = end_u_ - start_u_;
	int dy = end_v_ - start_v_;
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;

	double max_du = 0;
	double max_dv = 0;
	for(Basisfunction *b : basis_) {
		double du    = (*b)[0][order_u_] - (*b)[0][0];
		double dv    = (*b)[1][order_v_] - (*b)[1][0];
		max_du       = (max_du > du) ? max_du : du;
		max_dv       = (max_dv > dv) ? max_dv : dv;
	}

	double scaleSize = (max_du > max_dv) ? 1.0/max_du : 1.0/max_dv;
	scaleSize *= 20.0;
	
	// create the ellipse function
	out << "/ellipse {\n";
	out << "/endangle exch def\n";
	out << "/startangle exch def\n";
	out << "/yrad exch def\n";
	out << "/xrad exch def\n";
	out << "/y exch def\n";
	out << "/x exch def\n";
	out << "/savematrix matrix currentmatrix def\n";
	out << "x y translate\n";
	out << "xrad yrad scale\n";
	out << "0 0 1 startangle endangle arc\n";
	out << "savematrix setmatrix\n";
	out << "} def\n";

	// load the font to use
	out << "/Times-Roman findfont\n";
	out << "24 scalefont\n";
	out << "setfont\n";

	int i = -1;
	for(Basisfunction *b : basis_) {
		i++;
		double avg_u = 0;
		double avg_v = 0;
		double du    = (*b)[0][order_u_] - (*b)[0][0];
		double dv    = (*b)[1][order_v_] - (*b)[1][0];

		// move C^{-1} text on internal functions
		double textOffset = 15.0;
		int textX   = ((*b)[0][1] == (*b)[0][order_u_]) ? -2 : 1;
		int textY   = ((*b)[1][1] == (*b)[1][order_v_]) ? -2 : 1;
		// move text on edge functions
		if((*b)[0][1] == end_u_)
			textX = 1;
		else if((*b)[0][order_u_-1] == start_u_)
			textX = -2;
		if((*b)[1][1] == end_v_)
			textY = 1;
		else if((*b)[1][order_v_-1] == start_v_)
			textY = -2;

		for(int j=1; j<order_u_; j++)
			avg_u += (*b)[0][j];
		for(int j=1; j<order_v_; j++)
			avg_v += (*b)[1][j];
		avg_u /= (order_u_-1);
		avg_v /= (order_v_-1);

		bool doColor = false;
		if(colorBasis != NULL)
			for(uint j=0; j<colorBasis->size(); j++)
				if(i == (uint) colorBasis->at(j))
					doColor = true;

		if(doColor) {
			out << selected_basis_red   << " ";
			out << selected_basis_green << " ";
			out << selected_basis_blue  << " ";
			out << "setrgbcolor \n";
		} else {
			out << basis_red   << " ";
			out << basis_green << " ";
			out << basis_blue  << " ";
			out << "setrgbcolor \n";
		}
		if(drawAll || doColor) {
			out << "newpath\n";
			out << avg_u*scale << " " << avg_v*scale << " " << du*scaleSize << " " << dv*scaleSize << " 0 360 ellipse\n";
			out << "closepath fill\n";
			out << "0 setgray\n";
			out << avg_u*scale << " " << avg_v*scale << " " << du*scaleSize << " " << dv*scaleSize << " 0 360 ellipse\n";
			out << "closepath stroke\n";
		}

		out << "\n";
		out << "newpath\n";
		out << avg_u*scale + textX*textOffset << " " << avg_v*scale + textY*textOffset << " moveto\n";
		out << "(" << i << ") show\n";
		out << "\n";
	}
	if(close)
		out << "%%EOF\n";
}

void LRSplineSurface::printElements(std::ostream &out) const {
	for(uint i=0; i<element_.size(); i++) {
		if(i<100) out << " ";
		if(i<10)  out << " ";
		out << i << ": " << *element_[i] << std::endl;
	}
}

#undef MY_STUPID_FABS
#undef DOUBLE_TOL

} // end namespace LR

