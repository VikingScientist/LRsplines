#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/MeshRectangle.h"
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

#define TIME_LRSPLINE 1


typedef unsigned int uint;

//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;

namespace LR {

#define DOUBLE_TOL 1e-14
#define MY_STUPID_FABS(x) (((x)>0)?(x):-(x))


LRSplineVolume::LRSplineVolume() {
	rational_ = false;
	dim_      = 0;
	order_u_  = 0;
	order_v_  = 0;
	order_w_  = 0;
	start_u_  = 0;
	start_v_  = 0;
	start_w_  = 0;
	end_u_    = 0;
	end_v_    = 0;
	end_w_    = 0;
	meshrect_ = std::vector<MeshRectangle*>(0);
	element_  = std::vector<Element*>(0);
	maxTjoints_       = -1;
	doCloseGaps_      = true;
	maxAspectRatio_   = 2.0;
	doAspectRatioFix_ = false;
	refStrat_         = LR_SAFE;
	refKnotlineMult_  = 1;
	symmetry_         = 1;
}

LRSplineVolume::LRSplineVolume(Go::SplineVolume *vol) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	rational_ = vol->rational();
	dim_      = vol->dimension();
	order_u_  = vol->order(0);
	order_v_  = vol->order(1);
	order_w_  = vol->order(2);
	start_u_  = vol->startparam(0);
	start_v_  = vol->startparam(1);
	start_w_  = vol->startparam(2);
	end_u_    = vol->endparam(0);
	end_v_    = vol->endparam(1);
	end_w_    = vol->endparam(2);
	maxTjoints_       = -1;
	doCloseGaps_      = true;
	maxAspectRatio_   = 2.0;
	doAspectRatioFix_ = false;
	refStrat_         = LR_SAFE;
	refKnotlineMult_  = 1;
	symmetry_         = 1;

	int n1 = vol->numCoefs(0);
	int n2 = vol->numCoefs(1);
	int n3 = vol->numCoefs(2);
	std::vector<double>::iterator coef  = vol->coefs_begin();
	std::vector<double>::const_iterator knot_u = vol->basis(0).begin();
	std::vector<double>::const_iterator knot_v = vol->basis(1).begin();
	std::vector<double>::const_iterator knot_w = vol->basis(2).begin();
	for(int k=0; k<n3; k++)
		for(int j=0; j<n2; j++)
			for(int i=0; i<n1; i++)
				basis_.insert(new Basisfunction(&(*(knot_u+i)),
				                                &(*(knot_v+j)),
				                                &(*(knot_w+k)),
				                                &(*(coef+(k*n1*n2+j*n1+i)*(dim_+rational_))),
				                                dim_, order_u_, order_v_, order_w_) );
	int unique_u=0;
	int unique_v=0;
	int unique_w=0;
	for(int i=0; i<n1+order_u_; i++) {// const u, spanning v
		int mult = 1;
		while(i+1<n1+order_u_ && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],  knot_w[0],
		                                      knot_u[i], knot_v[n2], knot_w[n3], mult));
	}
	for(int i=0; i<n2+order_v_; i++) {// const v, spanning u
		int mult = 1;
		while(i+1<n2+order_v_ && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[i], knot_w[0],
		                                      knot_u[n1], knot_v[i], knot_w[n3], mult));
	}
	for(int i=0; i<n3+order_w_; i++) {
		int mult = 1;
		while(i+1<n3+order_w_ && knot_w[i]==knot_w[i+1]) {
			i++;
			mult++;
		}
		unique_w++;
		meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[0],  knot_w[i],
		                                      knot_u[n1], knot_v[n2], knot_w[i], mult));
	}
	for(int k=0; k<unique_w-1; k++) {
		for(int j=0; j<unique_v-1; j++) {
			for(int i=0; i<unique_u-1; i++) {
				double umin = meshrect_[                      i  ]->start_[0];
				double vmin = meshrect_[unique_u +            j  ]->start_[1];
				double wmin = meshrect_[unique_v + unique_u + k  ]->start_[2];
				double umax = meshrect_[                      i+1]->stop_[0];
				double vmax = meshrect_[unique_u +            j+1]->stop_[1];
				double wmax = meshrect_[unique_v + unique_u + k+1]->stop_[2];
				double min[] = {umin, vmin, wmin};
				double max[] = {umax, vmax, wmax};
				element_.push_back(new Element(3, min, max));
			}
		}
	}
	for(Basisfunction* b : basis_)
		updateSupport(b);
}


LRSplineVolume::LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, double *knot_u, double *knot_v, double *knot_w, double *coef, int dim, bool rational) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	rational_ = rational;
	dim_      = dim;
	order_u_  = order_u;
	order_v_  = order_v;
	order_w_  = order_w;
	start_u_  = knot_u[0];
	start_v_  = knot_v[0];
	start_w_  = knot_w[0];
	end_u_    = knot_u[n1];
	end_v_    = knot_v[n2];
	end_w_    = knot_w[n3];
	maxTjoints_       = -1;
	doCloseGaps_      = true;
	maxAspectRatio_   = 2.0;
	doAspectRatioFix_ = false;
	refStrat_         = LR_SAFE;
	refKnotlineMult_  = 1;
	symmetry_         = 1;

	for(int k=0; k<n3; k++)
		for(int j=0; j<n2; j++)
			for(int i=0; i<n1; i++)
				basis_.insert(new Basisfunction(knot_u+i,
				                                knot_v+j,
				                                knot_w+k,
				                                coef+(k*n1*n2+j*n1+i)*(dim_+rational_),
				                                dim_, order_u_, order_v_, order_w_) );
	int unique_u=0;
	int unique_v=0;
	int unique_w=0;
	for(int i=0; i<n1+order_u_; i++) {// const u, spanning v
		int mult = 1;
		while(i+1<n1+order_u_ && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],  knot_w[0],
		                                      knot_u[i], knot_v[n2], knot_w[n3], mult));
	}
	for(int i=0; i<n2+order_v_; i++) {// const v, spanning u
		int mult = 1;
		while(i+1<n2+order_v_ && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[i], knot_w[0],
		                                      knot_u[n1], knot_v[i], knot_w[n3], mult));
	}
	for(int i=0; i<n3+order_w_; i++) {
		int mult = 1;
		while(i+1<n3+order_w_ && knot_w[i]==knot_w[i+1]) {
			i++;
			mult++;
		}
		unique_w++;
		meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[0],  knot_w[i],
		                                      knot_u[n1], knot_v[n2], knot_w[i], mult));
	}
	for(int k=0; k<unique_w-1; k++) {
		for(int j=0; j<unique_v-1; j++) {
			for(int i=0; i<unique_u-1; i++) {
				double umin = meshrect_[                      i  ]->start_[0];
				double vmin = meshrect_[unique_u +            j  ]->start_[1];
				double wmin = meshrect_[unique_v + unique_u + k  ]->start_[2];
				double umax = meshrect_[                      i+1]->stop_[0];
				double vmax = meshrect_[unique_u +            j+1]->stop_[1];
				double wmax = meshrect_[unique_v + unique_u + k+1]->stop_[2];
				double min[] = {umin, vmin, wmin};
				double max[] = {umax, vmax, wmax};
				element_.push_back(new Element(3, min, max));
			}
		}
	}
	for(Basisfunction* b : basis_)
		updateSupport(b);

}

LRSplineVolume::~LRSplineVolume() {
	for(Basisfunction* b : basis_)
		delete b;
	for(uint i=0; i<meshrect_.size(); i++)
		delete meshrect_[i];
	for(uint i=0; i<element_.size(); i++)
		delete element_[i];
}


LRSplineVolume* LRSplineVolume::copy()
{
	/*
	LRSplineVolume *returnvalue = new LR::LRSplineVolume();
	
	for(int i = 0; i< this->nBasisFunctions();i++)
		{returnvalue -> basis_.push_back(this->basis_[i]->copy());}
	
	for(int i = 0; i < this->nElements();i++)
		{returnvalue -> element_.push_back(this->element_[i]->copy());}
	
	for(int i = 0; i < this->nMeshRectangles();i++)
		{returnvalue -> meshrect_.push_back(this-> meshrect_[i]->copy());}
	
	returnvalue->rational_         = this->rational_;
	returnvalue->dim_              = this->dim_;
	returnvalue->order_u_          = this->order_u_;
	returnvalue->order_v_          = this->order_v_;
	returnvalue->order_w_          = this->order_w_;
	returnvalue->start_u_          = this->start_u_;
	returnvalue->start_v_          = this->start_v_;
	returnvalue->start_w_          = this->start_w_;
	returnvalue->end_u_            = this->end_u_;
	returnvalue->end_v_            = this->end_v_;
	returnvalue->end_w_            = this->end_w_;
	returnvalue->maxTjoints_       = this->maxTjoints_;
	returnvalue->doCloseGaps_      = this->doCloseGaps_;
	returnvalue->doAspectRatioFix_ = this->doAspectRatioFix_;
	returnvalue->maxAspectRatio_   = this->maxAspectRatio_;

	
	for(int i = 0; i< this->nBasisFunctions();i++)
		{returnvalue -> updateSupport(returnvalue->basis_[i]);}
	
	return returnvalue;
	*/
	return NULL;
}

LRSplineVolume& LRSplineVolume::operator=( LRSplineVolume &copythis) {

	this->basis_     = copythis.basis_;
	this->rational_  = copythis.rational_;
	
	this->meshrect_  = copythis.meshrect_;
	this->element_   = copythis.element_;
	this->dim_       = copythis.dim_;
	this->order_u_   = copythis.order_u_;
	this->order_v_   = copythis.order_v_;
	this->order_w_   = copythis.order_w_;
	this->start_u_   = copythis.start_u_;
	this->start_v_   = copythis.start_v_;
	this->start_w_   = copythis.start_w_;
	this->end_u_     = copythis.end_u_;
	this->end_v_     = copythis.end_v_;
	this->end_w_     = copythis.end_w_;

	return *this;
}

void LRSplineVolume::point(Go::Point &pt, double u, double v, double w, int iEl, bool u_from_right, bool v_from_right, bool w_from_right) const {
	Go::Point cp;
	double basis_ev;

	pt.resize(dim_);
	pt.setValue(0.0);

	for(Basisfunction* b : element_[iEl]->support()) {
		b->getControlPoint(cp);
		
		basis_ev = b->evaluate(u,v,w, u_from_right, v_from_right, w_from_right);
		pt += basis_ev*cp;
	}
	
}

void LRSplineVolume::point(Go::Point &pt, double u, double v, double w, int iEl) const {
	Go::Point cp;
	double basis_ev;

	pt.resize(dim_);
	pt.setValue(0.0);
	if(iEl == -1)
		iEl = getElementContaining(u,v,w);
	
	for(Basisfunction* b : element_[iEl]->support()) {
		b->getControlPoint(cp);
		
		basis_ev = b->evaluate(u,v,w,  u!=end_u_, v!=end_v_, w!=end_w_);
		pt += basis_ev*cp;
	}
	
}

void LRSplineVolume::point(std::vector<Go::Point> &pts, double u, double v, double w, int derivs, int iEl) const {
#ifdef TIME_LRSPLINE
	PROFILE("Point()");
#endif
	Go::Point cp;
	std::vector<double> basis_ev;

	// clear and resize output array (optimization may consider this an outside task)
	pts.resize((derivs+1)*(derivs+2)*(2*derivs+6)/12);
	for(uint i=0; i<pts.size(); i++) {
		pts[i].resize(dim_);
		pts[i].setValue(0.0);
	}

	if(iEl == -1)
		iEl = getElementContaining(u,v,w);

	for(Basisfunction* b : element_[iEl]->support() ) {
		b->getControlPoint(cp);
		b->evaluate(basis_ev, u,v,w, derivs, u!=end_u_, v!=end_v_, w!=end_w_);
		for(uint j=0; j<pts.size(); j++)
			pts[j] += basis_ev[j]*cp;
	}
}

void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w, Go::BasisDerivsSf2 & result, int iEl ) const {
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
		(*it)->evaluate(values, param_u, param_v, 2, param_u!=end_u_, param_v!=end_v_);
	
		result.basisValues[i]    = values[0];
		result.basisDerivs_u[i]  = values[1];
		result.basisDerivs_v[i]  = values[2];
		result.basisDerivs_uu[i] = values[3];
		result.basisDerivs_uv[i] = values[4];
		result.basisDerivs_vv[i] = values[5];
	}
}

void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w, Go::BasisDerivsSf & result, int iEl ) const {
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


void LRSplineVolume::computeBasis(double param_u, double param_v, double param_w, Go::BasisPtsSf & result, int iEl ) const {
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.preparePts(param_u, param_v, 0, -1, nPts);
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i)
		result.basisValues[i] = (*it)->evaluate(param_u, param_v, param_u!=end_u_, param_v!=end_v_, param_w!=end_w_);
}

void LRSplineVolume::computeBasis (double param_u,
                                   double param_v,
                                   double param_w,
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
	    (*it)->evaluate(result[i], param_u, param_v, derivs, param_u!=end_u_, param_v!=end_v_, param_w!=end_w_);
}

int LRSplineVolume::getElementContaining(double u, double v, double w) const {
	for(uint i=0; i<element_.size(); i++)
		if(element_[i]->getParmin(0) <= u && element_[i]->getParmin(1) <= v && element_[i]->getParmin(2) <= w) 
			if((u < element_[i]->getParmax(0) || (u == end_u_ && u <= element_[i]->getParmax(0))) && 
			   (v < element_[i]->getParmax(1) || (v == end_v_ && v <= element_[i]->getParmax(1))) && 
			   (w < element_[i]->getParmax(1) || (w == end_w_ && w <= element_[i]->getParmax(2))))
				return i;
		
	return -1;
}

/*
void LRSplineVolume::getMinspanLines(int iEl, std::vector<MeshRectangle*>& lines) {
	Element *e = element_[iEl];
	std::vector<Basisfunction*>::iterator it;
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
	for(it=e->supportBegin(); it<e->supportEnd(); it++) {
		double lowu  = (**it).umin();
		double highu = (**it).umax();
		double lowv  = (**it).vmin();
		double highv = (**it).vmax();
		double du = highu - lowu;
		double dv = highv - lowv;
		int startI=0;
		int stopI=0;
		int startJ=0;
		int stopJ=0;
		while((**it).knot_u_[startI] <= e->umin())
			startI++;
		while((**it).knot_u_[stopI]  <  e->umax())
			stopI++;
		while((**it).knot_v_[startJ] <= e->vmin())
			startJ++;
		while((**it).knot_v_[stopJ]  <  e->vmax())
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
		lines.push_back(new MeshRectangle(true, (e->vmin() + e->vmax())/2.0, umin, umax, 1));
		
	if(!only_insert_span_u_line) 
		lines.push_back(new MeshRectangle(false, (e->umin() + e->umax())/2.0, vmin, vmax, 1));

}

void LRSplineVolume::getFullspanLines(int iEl, std::vector<MeshRectangle*>& lines) {
	std::vector<Basisfunction*>::iterator it;
	Element *e = element_[iEl];
	double umin = e->umin();
	double umax = e->umax();
	double vmin = e->vmin();
	double vmax = e->vmax();
	bool   only_insert_span_u_line = (vmax-vmin) >= maxAspectRatio_*(umax-umin);
	bool   only_insert_span_v_line = (umax-umin) >= maxAspectRatio_*(vmax-vmin);
	// loop over all supported B-splines and make sure that everyone is covered by meshrect
	for(it=e->supportBegin(); it<e->supportEnd(); it++) {
		umin = (umin > (**it).umin()) ? (**it).umin() : umin;
		umax = (umax < (**it).umax()) ? (**it).umax() : umax;
		vmin = (vmin > (**it).vmin()) ? (**it).vmin() : vmin;
		vmax = (vmax < (**it).vmax()) ? (**it).vmax() : vmax;
	}
	if(!only_insert_span_v_line) 
		lines.push_back(new MeshRectangle(true, (e->vmin() + e->vmax())/2.0, umin, umax, 1));
		
	if(!only_insert_span_u_line) 
		lines.push_back(new MeshRectangle(false, (e->umin() + e->umax())/2.0, vmin, vmax, 1));
}

void LRSplineVolume::getStructMeshLines(int iBasis, std::vector<MeshRectangle*>& lines) {
	Basisfunction *b = basis_[iBasis];
	double umin = b->umin();
	double umax = b->umax();
	double vmin = b->vmin();
	double vmax = b->vmax();

	// find the largest knotspan in this function
	double max_du = 0;
	double max_dv = 0;
	for(int j=0; j<order_u_; j++) {
		double du = b->knot_u_[j+1]-b->knot_u_[j];
		bool isZeroSpan =  MY_STUPID_FABS(du) < DOUBLE_TOL ;
		max_du = (isZeroSpan || max_du>du) ? max_du : du;
	}
	for(int j=0; j<order_v_; j++) {
		double dv = b->knot_v_[j+1]-b->knot_v_[j];
		bool isZeroSpan =  MY_STUPID_FABS(dv) < DOUBLE_TOL ;
		max_dv = (isZeroSpan || max_dv>dv) ? max_dv : dv;
	}

	// to keep as "square" basis function as possible, only insert
	// into the largest knot spans
	for(int j=0; j<order_u_; j++) {
		double du = b->knot_u_[j+1]-b->knot_u_[j];
		if( MY_STUPID_FABS(du-max_du) < DOUBLE_TOL )
			lines.push_back(new MeshRectangle(false, (b->knot_u_[j] + b->knot_u_[j+1])/2.0, vmin, vmax,1));
	}
	for(int j=0; j<order_v_; j++) {
		double dv = b->knot_v_[j+1]-b->knot_v_[j];
		if( MY_STUPID_FABS(dv-max_dv) < DOUBLE_TOL )
			lines.push_back(new MeshRectangle(true, (b->knot_v_[j] + b->knot_v_[j+1])/2.0, umin, umax,1));
	}
}
*/

#if 0
void LRSplineVolume::refineBasisFunction(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineBasisFunction(tmp);
}

void LRSplineVolume::refineBasisFunction(const std::vector<int> &indices) {
	std::vector<MeshRectangle*> newLines;

	/* first retrieve all meshrects needed */
	for(uint i=0; i<indices.size(); i++)
		getStructMeshLines(indices[i],newLines);

	/* Do the actual refinement */
	for(uint i=0; i<newLines.size(); i++) {
		MeshRectangle *m = newLines[i];
		insert_line(!m->is_spanning_u(), m->const_par_, m->start_, m->stop_, refKnotlineMult_);
	}

	/* do a posteriori fixes to ensure a proper mesh */
	aPosterioriFixes();

	/* exit cleanly be deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++) 
		delete newLines[i];
}

void LRSplineVolume::refineElement(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineElement(tmp);
}

void LRSplineVolume::refineElement(const std::vector<int> &indices) {
	std::vector<MeshRectangle*> newLines;

	/* first retrieve all meshrects needed */
	for(uint i=0; i<indices.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanLines(indices[i],newLines);
		else
			getFullspanLines(indices[i],newLines);
	}

	/* Do the actual refinement */
	for(uint i=0; i<newLines.size(); i++) {
		MeshRectangle *m = newLines[i];
		insert_line(!m->is_spanning_u(), m->const_par_, m->start_, m->stop_, refKnotlineMult_);
	}

	/* do a posteriori fixes to ensure a proper mesh */
	aPosterioriFixes();

	/* exit cleanly be deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++) 
		delete newLines[i];
}

void LRSplineVolume::refineByDimensionIncrease(const std::vector<double> &errPerElement, double beta) {
	Basisfunction *b;
	Element       *e;
	/* accumulate the error & index - vector */
	std::vector<IndexDouble> errors;
	if(refStrat_ == LR_ISOTROPIC_FUNC) { // error per-function
		for(uint i=0; i<basis_.size(); i++) {
			b = basis_[i];
			errors.push_back(IndexDouble(0.0, i));
			for(int j=0; j<b->nSupportedElements(); j++) {
				e = b->support_[j];
				errors[i].first += errPerElement[e->getId()];
			}
		}
	} else {
		for(uint i=0; i<element_.size(); i++) 
			errors.push_back(IndexDouble(errPerElement[i], i));
	}

	/* sort errors */
	std::sort(errors.begin(), errors.end(), std::greater<IndexDouble>());

	/* first retrieve all possible meshrects needed */
	std::vector<std::vector<MeshRectangle*> > newLines(errors.size(), std::vector<MeshRectangle*>(0));
	for(uint i=0; i<errors.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanLines(errors[i].second, newLines[i]);
		else if(refStrat_ == LR_SAFE) 
			getFullspanLines(errors[i].second, newLines[i]);
		else if(refStrat_ == LR_ISOTROPIC_FUNC) 
			getStructMeshLines(errors[i].second, newLines[i]);
		// note that this is an excessive loop as it computes the meshrects for ALL elements,
		// but we're only going to use a small part of this.
	}

	/* Do the actual refinement */
	uint target_n_functions = ceil(basis_.size()*(1+beta));
	int i=0;
	while( basis_.size() < target_n_functions ) {
		for(uint j=0; j<newLines[i].size(); j++) {
			MeshRectangle *m = newLines[i][j];
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

void LRSplineVolume::aPosterioriFixes()  {
	std::vector<MeshRectangle*> *newLines = NULL;
	uint nFunc;
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


void LRSplineVolume::closeGaps(std::vector<MeshRectangle*>* newLines) {
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
		for(uint j=0; j<meshrect_.size(); j++) {
			MeshRectangle *m = meshrect_[j];
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
	MeshRectangle* m;
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

void LRSplineVolume::enforceMaxTjoints(std::vector<MeshRectangle*> *newLines) {
	bool someFix = true;
	while(someFix) {
		someFix = false;
		for(uint i=0; i<element_.size(); i++) {
			double umin = element_[i]->umin();
			double umax = element_[i]->umax();
			double vmin = element_[i]->vmin();
			double vmax = element_[i]->vmax();
			std::vector<double> left, right, top, bottom;
			for(uint j=0; j<meshrect_.size(); j++) {
				MeshRectangle *m = meshrect_[j];
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
			MeshRectangle *m;
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

void LRSplineVolume::enforceMaxAspectRatio(std::vector<MeshRectangle*>* newLines) {
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
				std::vector<MeshRectangle*> splitLines; // should always contain exactly one meshrect on function return
				if(refStrat_ == LR_MINSPAN) 
					getMinspanLines(i, splitLines);
				else
					getFullspanLines(i, splitLines);

				
				MeshRectangle *m, *msplit;
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

MeshRectangle* LRSplineVolume::insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity) {
	return insert_line(true, u, start_v, stop_v, multiplicity);
}
#endif 

MeshRectangle* LRSplineVolume::insert_line(MeshRectangle *newRect) {
	if(newRect->start_[0] < start_u_ ||
	   newRect->start_[1] < start_v_ ||
	   newRect->start_[2] < start_w_ ||
	   newRect->stop_[0]  > end_u_  ||
	   newRect->stop_[1]  > end_v_  ||
	   newRect->stop_[2]  > end_w_ ) {
		std::cerr << "Error: inserting meshrctangle " << *newRect << " outside parametric domain";
		std::cerr << "(" << start_u_ << ", " << start_v_ << ", " << start_w_ << ") x ";
		std::cerr << "(" <<   end_u_ << ", " <<   end_v_ << ", " <<   end_w_ << ") ";
		return NULL;
	}


	{ // check if the line is an extension or a merging of existing lines
#ifdef TIME_LRSPLINE
	PROFILE("meshrectangle verification");
#endif
	for(uint i=0; i<meshrect_.size(); i++) {
		// if newRect overlaps any existing ones (may be multiple existing ones)
		// let newRect be the entire length of all merged and delete the unused ones
		if(meshrect_[i]->overlaps(newRect)) {
			if(meshrect_[i]->equals(newRect)) {
				if(meshrect_[i]->multiplicity_ >= newRect->multiplicity_) {
					return meshrect_[i];
				} else {
					// keeping newrect, getting rid of the old rect
					delete meshrect_[i];
					meshrect_.erase(meshrect_.begin() + i);
					i--;
				}
			} else {
				std::cerr << "Haven't fixed overlapping, nonequal meshrectangles yet\n";
				std::cerr << "Just works with brand new shiny rectangles\n";
				exit(54023993);
			}
		}
	}
	} // end meshrectangle verification timer

	HashSet<Basisfunction*> newFuncStp1, newFuncStp2;
	HashSet<Basisfunction*> removeFunc;

	{ // STEP 1: test EVERY function against the NEW meshrect
#ifdef TIME_LRSPLINE
	PROFILE("STEP 1");
#endif
	for(Basisfunction* b : basis_) {
		if(newRect->splits(b)) {
			int nKnots;
			if( ((nKnots=newRect->nKnotsIn(b)) != newRect->multiplicity_) ) {
				removeFunc.insert(b);
				split( newRect->constDirection(), b, newRect->constParameter(), newRect->multiplicity_-nKnots, newFuncStp1 ); 
			}
		}
	}
	for(Basisfunction* b : removeFunc) {
		basis_.erase(b);
		delete b;
	}
	for(uint i=0; i<element_.size(); i++) {
		if(newRect->splits(element_[i]))
			element_.push_back(element_[i]->split(newRect->constDirection(), newRect->constParameter()) );
	}
	} // end step 1 timer

	{ // STEP 2: test every NEW function against ALL old meshrects
#ifdef TIME_LRSPLINE
	PROFILE("STEP 2");
#endif
	meshrect_.push_back(newRect);
	while(newFuncStp1.size() > 0) {
		Basisfunction *b = newFuncStp1.pop();
		bool splitMore = false;
		for(MeshRectangle *m : meshrect_) {
			if(m->splits(b)) {
				int nKnots = m->nKnotsIn(b);
				if( nKnots != m->multiplicity_ ) {
					splitMore = true;
					split( m->constDirection(), b, m->constParameter(), m->multiplicity_-nKnots, newFuncStp1);
					delete b;
					break;
				}
			}
		}
		if(!splitMore)
			basis_.insert(b);
	}
	} // end step 2 timer
	return newRect;
}

#if 0

MeshRectangle* LRSplineVolume::insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity) {
	return insert_line(false, v, start_u, stop_u, multiplicity);
}

#endif
 

int LRSplineVolume::split(int constDir, Basisfunction *b, double new_knot, int multiplicity, HashSet<Basisfunction*> &newFunctions) {
#ifdef TIME_LRSPLINE
	PROFILE("split()");
#endif

	// create the new functions b1 and b2
	Basisfunction *b1, *b2;
	std::vector<double> knot = b->getknots(constDir);
	int     p                = b->getOrder(constDir);
	int     insert_index     = 0;
	if(new_knot < knot[0] || knot[p] < new_knot)
		return 0;
	while(knot[insert_index] < new_knot)
		insert_index++;
	double alpha1 = (insert_index == p)  ? 1.0 : (new_knot-knot[0])/(knot[p-1]-knot[0]);
	double alpha2 = (insert_index == 1 ) ? 1.0 : (knot[p]-new_knot)/(knot[p]-knot[1]);
	double newKnot[p+2];
	std::copy(knot.begin(), knot.begin()+p+1, newKnot+1);
	newKnot[0] = new_knot;
	std::sort(newKnot, newKnot + p+2);
	if(constDir == 0) {
		b1 = new Basisfunction(newKnot  ,  (*b)[1].begin(),  (*b)[2].begin(), b->cp(), b->dim(), order_u_, order_v_, order_w_, b->w()*alpha1);
		b2 = new Basisfunction(newKnot+1,  (*b)[1].begin(),  (*b)[2].begin(), b->cp(), b->dim(), order_u_, order_v_, order_w_, b->w()*alpha2);
	} else if(constDir == 1) {
		b1 = new Basisfunction((*b)[0].begin(), newKnot   ,  (*b)[2].begin(), b->cp(), b->dim(), order_u_, order_v_, order_w_, b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), newKnot+1 ,  (*b)[2].begin(), b->cp(), b->dim(), order_u_, order_v_, order_w_, b->w()*alpha2);
	} else { // insert in w
		b1 = new Basisfunction((*b)[0].begin(), (*b)[1].begin(),  newKnot   , b->cp(), b->dim(), order_u_, order_v_, order_w_, b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), (*b)[1].begin(),  newKnot+1 , b->cp(), b->dim(), order_u_, order_v_, order_w_, b->w()*alpha2);
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
			bool recursive_split = (multiplicity > 1) && (*b1)[constDir].back() != new_knot;
			if(recursive_split) {
				split( constDir, b1, new_knot, multiplicity-1, newFunctions);
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
			bool recursive_split = (multiplicity > 1) && (*b2)[constDir][0] != new_knot;
			if(recursive_split) {
				split( constDir, b2, new_knot, multiplicity-1, newFunctions);
				delete b2;
			} else {
				newFunctions.insert(b2);
			}
		}
	}
}

void LRSplineVolume::getEdgeFunctions(std::vector<Basisfunction*> &edgeFunctions, parameterEdge edge, int depth) const {
	edgeFunctions.clear();
	for(Basisfunction *b : basis_) {
		bool ok = true;
		if( edge & WEST )
			if((*b)[0][order_u_-depth] != start_u_) 
				ok = false;
		if( edge & EAST )
			if((*b)[0][depth] != end_u_) 
				ok = false;
		if( edge & SOUTH )
			if((*b)[1][order_v_-depth] != start_v_) 
				ok = false;
		if( edge & NORTH )
			if((*b)[1][depth] != end_v_) 
				ok = false;
		if( edge & BOTTOM )
			if((*b)[2][order_w_-depth] != start_w_) 
				ok = false;
		if( edge & TOP )
			if((*b)[2][depth] != end_u_) 
				ok = false;

		if(ok)
			edgeFunctions.push_back(b);
	}
}

void LRSplineVolume::rebuildDimension(int dimvalue) {
	for(Basisfunction* b : basis_)
		b->setDimension(dimvalue);
	dim_ = dimvalue;
}

#if 0
void LRSplineVolume::getGlobalKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
	getGlobalUniqueKnotVector(knot_u, knot_v);

	// add in duplicates where apropriate
	for(uint i=0; i<knot_u.size(); i++) {
		for(uint j=0; j<meshrect_.size(); j++) {
			if(!meshrect_[j]->is_spanning_u() && meshrect_[j]->const_par_==knot_u[i]) {
				for(int m=1; m<meshrect_[j]->multiplicity_; m++) {
					knot_u.insert(knot_u.begin()+i, knot_u[i]);
					i++;
				}
				break;
			}
		}
	}
	for(uint i=0; i<knot_v.size(); i++) {
		for(uint j=0; j<meshrect_.size(); j++) {
			if(meshrect_[j]->is_spanning_u() && meshrect_[j]->const_par_==knot_v[i]) {
				for(int m=1; m<meshrect_[j]->multiplicity_; m++) {
					knot_v.insert(knot_v.begin()+i, knot_v[i]);
					i++;
				}
				break;
			}
		}
	}
}

void LRSplineVolume::getGlobalUniqueKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
	knot_u.clear();
	knot_v.clear();
	// create a huge list of all line instances
	for(uint i=0; i<meshrect_.size(); i++) {
		if(meshrect_[i]->is_spanning_u())
			knot_v.push_back(meshrect_[i]->const_par_);
		else
			knot_u.push_back(meshrect_[i]->const_par_);
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
#endif

void LRSplineVolume::updateSupport(Basisfunction *f,
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

void LRSplineVolume::updateSupport(Basisfunction *f) {
	updateSupport(f, element_.begin(), element_.end());
}

void LRSplineVolume::generateIDs() const {
	uint i=0;
	for(Basisfunction *b : basis_)
		b->setId(i++);
	for(i=0; i<element_.size(); i++) 
		element_[i]->setId(i);
}

bool LRSplineVolume::isLinearIndepByOverloading(bool verbose) {
	std::vector<Basisfunction*>           overloaded;
	std::vector<Element*>::iterator       eit;
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
			for(eit=overloaded[i]->supportedElementBegin(); eit<overloaded[i]->supportedElementEnd(); eit++)
				(*eit)->incrementOverloadCount();
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
		}


		iterationCount++;
	} while((uint) lastOverloadCount != overloaded.size());

	return overloaded.size() == 0;
}

void LRSplineVolume::getBezierElement(int iEl, std::vector<double> &controlPoints) const {
	controlPoints.clear();
	controlPoints.resize(order_u_*order_v_*order_w_*dim_, 0);
	Element *el = element_[iEl];
	for(Basisfunction* b : el->support()) {
		std::vector<double> knotU((*b)[0].begin(), (*b)[0].end() );
		std::vector<double> knotV((*b)[1].begin(), (*b)[1].end() );
		std::vector<double> knotW((*b)[2].begin(), (*b)[2].end() );
		int startU=-1;
		int startV=-1;
		int startW=-1;

		std::vector<double> rowU(1,1), rowV(1,1), rowW(1,1);

		double min = el->umin();
		double max = el->umax();
		while(knotU[++startU] < min);
		// if(knotU[startU+order_u_-1] == min) startU++;
		while(true) {
			int p    = order_u_-1;
			int newI = -1;
			double z;
			if(       knotU.size() < startU+order_u_   || knotU[startU+  order_u_-1] != min) {
				z    = min;
				newI = startU;
			} else if(knotU.size() < startU+2*order_u_ || knotU[startU+2*order_u_-1] != max ) {
				z    = max;
				newI = startU + order_u_;
			} else {
				break;
			}
		      
			std::vector<double> newRowU(rowU.size()+1, 0);
			for(uint k=0; k<rowU.size(); k++) {
				#define U(x) ( knotU[x+k] )
				if(z < U(0) || z > U(p+1)) {
					newRowU[k] = rowU[k];
					continue;
				}
				double alpha1 = (U(p) <=  z  ) ? 1 : double(   z    - U(0)) / (  U(p)  - U(0));
				double alpha2 = (z    <= U(1)) ? 1 : double( U(p+1) - z   ) / ( U(p+1) - U(1));
				newRowU[k]   += rowU[k]*alpha1;
				newRowU[k+1] += rowU[k]*alpha2;
				#undef U
			}
			knotU.insert(knotU.begin()+newI, z);
			rowU = newRowU;
		}
	
		min = el->vmin();
		max = el->vmax();
		while(knotV[++startV] < min);
		// if(knotV[startV+order_v_-1] == min) startV++;
		while(true) {
			int p    = order_v_-1;
			int newI = -1;
			double z;
			if(       knotV.size() < startV+order_v_   || knotV[startV+  order_v_-1] != min) {
				z    = min;
				newI = startV;
			} else if(knotV.size() < startV+2*order_v_ || knotV[startV+2*order_v_-1] != max ) {
				z = max;
				newI = startV + order_v_;
			} else {
				break;
			}
		      
			std::vector<double> newRowV(rowV.size()+1, 0);
			for(uint k=0; k<rowV.size(); k++) {
				#define V(x) ( knotV[x+k] )
				if(z < V(0) || z > V(p+1)) {
					newRowV[k] = rowV[k];
					continue;
				}
				double alpha1 = (V(p) <=  z  ) ? 1 : double(   z    - V(0)) / (  V(p)  - V(0));
				double alpha2 = (z    <= V(1)) ? 1 : double( V(p+1) - z   ) / ( V(p+1) - V(1));
				newRowV[k]   += rowV[k]*alpha1;
				newRowV[k+1] += rowV[k]*alpha2;
				#undef V
			}
			knotV.insert(knotV.begin()+newI, z);
			rowV = newRowV;
		}
	
		min = el->getParmin(2);
		max = el->getParmax(2);
		while(knotW[++startW] < min);
		// if(knotW[startW+order_w_-1] == min) startW++;
		while(true) {
			int p    = order_w_-1;
			int newI = -1;
			double z;
			if(       knotW.size() < startW+order_w_   || knotW[startW+  order_w_-1] != min) {
				z    = min;
				newI = startW;
			} else if(knotW.size() < startW+2*order_w_ || knotW[startW+2*order_w_-1] != max ) {
				z = max;
				newI = startW + order_w_;
			} else {
				break;
			}
		      
			std::vector<double> newRowW(rowW.size()+1, 0);
			for(uint k=0; k<rowW.size(); k++) {
				#define W(x) ( knotW[x+k] )
				if(z < W(0) || z > W(p+1)) {
					newRowW[k] = rowW[k];
					continue;
				}
				double alpha1 = (W(p) <=  z  ) ? 1 : double(   z    - W(0)) / (  W(p)  - W(0));
				double alpha2 = (z    <= W(1)) ? 1 : double( W(p+1) - z   ) / ( W(p+1) - W(1));
				newRowW[k]   += rowW[k]*alpha1;
				newRowW[k+1] += rowW[k]*alpha2;
				#undef W
			}
			knotW.insert(knotW.begin()+newI, z);
			rowW = newRowW;
		}
		
		int ip = 0;
		for(int w=startW; w<startW+order_w_; w++)
			for(int v=startV; v<startV+order_v_; v++)
				for(int u=startU; u<startU+order_u_; u++)
					for(int d=0; d<dim_; d++)
						controlPoints[ip++] += b->cp()[d]*rowU[u]*rowV[v]*rowW[w]*b->w();

	}
}

#if 0
bool LRSplineVolume::isLinearIndepByMappingMatrix(bool verbose) const {
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

	for (int i = 0; i < nmb_bas; ++i) {
		int startU, startV;
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_u_+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_v_+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_u_; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_v_; j++) {
			if(knots_v[startV] == basis_[i]->knot_v_[j]) startV--;
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

void LRSplineVolume::getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const {
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

	for (int i = 0; i < nmb_bas; ++i) {
		int startU, startV;
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_u_+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_v_+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_u_; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_v_; j++) {
			if(knots_v[startV] == basis_[i]->knot_v_[j]) startV--;
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
				Ct[(startV+i2)*n1 + (startU+i1)][i] = rowV[i2]*rowU[i1];
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

bool LRSplineVolume::isLinearIndepByFloatingPointMappingMatrix(bool verbose) const {
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

	for (int i = 0; i < nmb_bas; ++i) {
		int startU, startV;
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_u_+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_v_+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_u_; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_v_; j++) {
			if(knots_v[startV] == basis_[i]->knot_v_[j]) startV--;
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

double LRSplineVolume::makeIntegerKnots() {
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

	// scale all meshrect values according to this
	MeshRectangle *m;
	for(uint i=0; i<meshrect_.size(); i++) {
		m = meshrect_[i];
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
	Basisfunction *b;
	for(Basisfunction *b : basis_) {
		for(int j=0; j<order_u_+1; j++)
			b->knot_u_[j] = floor(b->knot_u_[j]/scale + 0.5);
		for(int j=0; j<order_v_+1; j++)
			b->knot_v_[j] = floor(b->knot_v_[j]/scale + 0.5);
	}
	
	// scale all LRSplineVolume values
	start_u_ = floor(start_u_/scale + 0.5);
	start_v_ = floor(start_v_/scale + 0.5);
	end_u_   = floor(end_u_  /scale + 0.5);
	end_v_   = floor(end_v_  /scale + 0.5);

	return scale;
}
#endif

void LRSplineVolume::getDiagonalElements(std::vector<int> &result) const  {
	result.clear();
	for(uint i=0; i<element_.size(); i++) 
		if(element_[i]->umin() == element_[i]->vmin() && element_[i]->umin() == element_[i]->getParmin(2))
			result.push_back(i);
}

void LRSplineVolume::getDiagonalBasisfunctions(std::vector<Basisfunction*> &result) const  {
	result.clear();
	for(Basisfunction *b : basis_) {
		bool isDiag = true;
		for(int j=0; j<order_u_+1; j++)
			if(b->getknots(0)[j] != b->getknots(1)[j] || b->getknots(0)[j] != b->getknots(2)[j])
				isDiag = false;
		if(isDiag)
			result.push_back(b);
	}
}


void LRSplineVolume::read(std::istream &is) {
	start_u_ =  DBL_MAX;
	end_u_   = -DBL_MAX;
	start_v_ =  DBL_MAX;
	end_v_   = -DBL_MAX;
	start_w_ =  DBL_MAX;
	end_w_   = -DBL_MAX;

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
	int nBasis, nElements, nMeshRectangles;
	is >> order_u_;    ws(is);
	is >> order_v_;    ws(is);
	is >> order_w_;    ws(is);
	is >> nBasis;      ws(is);
	is >> nMeshRectangles;  ws(is);
	is >> nElements;   ws(is);
	is >> dim_;        ws(is);
	is >> rational_;   ws(is);
	
	meshrect_.resize(nMeshRectangles);
	element_.resize(nElements);

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

	for(int i=0; i<nMeshRectangles; i++) {
		meshrect_[i] = new MeshRectangle();
		meshrect_[i]->read(is);
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
		start_u_ = (element_[i]->getParmin(0) < start_u_) ? element_[i]->getParmin(0) : start_u_;
		end_u_   = (element_[i]->getParmax(0) > end_u_  ) ? element_[i]->getParmax(0) : end_u_  ;
		start_v_ = (element_[i]->getParmin(1) < start_v_) ? element_[i]->getParmin(1) : start_v_;
		end_v_   = (element_[i]->getParmax(1) > end_v_  ) ? element_[i]->getParmax(1) : end_v_  ;
		start_w_ = (element_[i]->getParmin(2) < start_w_) ? element_[i]->getParmin(2) : start_w_;
		end_w_   = (element_[i]->getParmax(2) > end_w_  ) ? element_[i]->getParmax(2) : end_w_  ;
	}
}

void LRSplineVolume::write(std::ostream &os) const {
	generateIDs();
	os << std::setprecision(16);
	os << "# LRSPLINE\n";
	os << "#\tp1\tp2\tp3\tNbasis\tNline\tNel\tdim\trat\n\t";
	os << order_u_ << "\t";
	os << order_v_ << "\t";
	os << order_w_ << "\t";
	os << basis_.size() << "\t";
	os << meshrect_.size() << "\t";
	os << element_.size() << "\t";
	os << dim_ << "\t";
	os << rational_ << "\n";

	os << "# Basis functions:\n";
	for(Basisfunction* b : basis_) 
		os << *b << std::endl;
	os << "# Mesh rectangles:\n";
	for(MeshRectangle* m : meshrect_)  
		os << *m << std::endl;
	os << "# Elements:\n";
	for(Element* e : element_)
		os << *e << std::endl;
}

void LRSplineVolume::printElements(std::ostream &out) const {
	for(uint i=0; i<element_.size(); i++) {
		if(i<100) out << " ";
		if(i<10)  out << " ";
		out << i << ": " << *element_[i] << std::endl;
	}
}

#undef MY_STUPID_FABS
#undef DOUBLE_TOL

} // end namespace LR

