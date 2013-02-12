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


LRSplineVolume::LRSplineVolume() {
	order_.resize(3);
	start_.resize(3);
	end_.resize(3);
	rational_ = false;
	dim_      = 0;
	order_[0]  = 0;
	order_[1]  = 0;
	order_[2]  = 0;
	start_[0]  = 0;
	start_[1]  = 0;
	start_[2]  = 0;
	end_[0]    = 0;
	end_[1]    = 0;
	end_[2]    = 0;
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
	order_.resize(3);
	start_.resize(3);
	end_.resize(3);
	rational_ = vol->rational();
	dim_      = vol->dimension();
	order_[0]  = vol->order(0);
	order_[1]  = vol->order(1);
	order_[2]  = vol->order(2);
	start_[0]  = vol->startparam(0);
	start_[1]  = vol->startparam(1);
	start_[2]  = vol->startparam(2);
	end_[0]    = vol->endparam(0);
	end_[1]    = vol->endparam(1);
	end_[2]    = vol->endparam(2);
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
				                                dim_, order_[0], order_[1], order_[2]) );
	int unique_u=0;
	int unique_v=0;
	int unique_w=0;
	for(int i=0; i<n1+order_[0]; i++) {// const u, spanning v
		int mult = 1;
		while(i+1<n1+order_[0] && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],  knot_w[0],
		                                      knot_u[i], knot_v[n2], knot_w[n3], mult));
	}
	for(int i=0; i<n2+order_[1]; i++) {// const v, spanning u
		int mult = 1;
		while(i+1<n2+order_[1] && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[i], knot_w[0],
		                                      knot_u[n1], knot_v[i], knot_w[n3], mult));
	}
	for(int i=0; i<n3+order_[2]; i++) {
		int mult = 1;
		while(i+1<n3+order_[2] && knot_w[i]==knot_w[i+1]) {
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
	order_.resize(3);
	start_.resize(3);
	end_.resize(3);
	rational_ = rational;
	dim_      = dim;
	order_[0]  = order_u;
	order_[1]  = order_v;
	order_[2]  = order_w;
	start_[0]  = knot_u[0];
	start_[1]  = knot_v[0];
	start_[2]  = knot_w[0];
	end_[0]    = knot_u[n1];
	end_[1]    = knot_v[n2];
	end_[2]    = knot_w[n3];
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
				                                dim_, order_[0], order_[1], order_[2]) );
	int unique_u=0;
	int unique_v=0;
	int unique_w=0;
	for(int i=0; i<n1+order_[0]; i++) {// const u, spanning v
		int mult = 1;
		while(i+1<n1+order_[0] && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],  knot_w[0],
		                                      knot_u[i], knot_v[n2], knot_w[n3], mult));
	}
	for(int i=0; i<n2+order_[1]; i++) {// const v, spanning u
		int mult = 1;
		while(i+1<n2+order_[1] && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[i], knot_w[0],
		                                      knot_u[n1], knot_v[i], knot_w[n3], mult));
	}
	for(int i=0; i<n3+order_[2]; i++) {
		int mult = 1;
		while(i+1<n3+order_[2] && knot_w[i]==knot_w[i+1]) {
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


LRSplineVolume* LRSplineVolume::copy() const {
	generateIDs();

	std::vector<Basisfunction*> basisVector;

	// flat list to make it quicker to update pointers from Basisfunction to Element and back again
	LRSplineVolume *returnvalue = new LR::LRSplineVolume();
	
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
	
	for(MeshRectangle *m : meshrect_)
		returnvalue -> meshrect_.push_back(m->copy());
	
	returnvalue->rational_         = this->rational_;
	returnvalue->dim_              = this->dim_;
	returnvalue->order_[0]          = this->order_[0];
	returnvalue->order_[1]          = this->order_[1];
	returnvalue->order_[2]          = this->order_[2];
	returnvalue->start_[0]          = this->start_[0];
	returnvalue->start_[1]          = this->start_[1];
	returnvalue->start_[2]          = this->start_[2];
	returnvalue->end_[0]            = this->end_[0];
	returnvalue->end_[1]            = this->end_[1];
	returnvalue->end_[2]            = this->end_[2];
	returnvalue->maxTjoints_       = this->maxTjoints_;
	returnvalue->doCloseGaps_      = this->doCloseGaps_;
	returnvalue->doAspectRatioFix_ = this->doAspectRatioFix_;
	returnvalue->maxAspectRatio_   = this->maxAspectRatio_;
	
	return returnvalue;
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
		
		basis_ev = b->evaluate(u,v,w,  u!=end_[0], v!=end_[1], w!=end_[2]);
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
		b->evaluate(basis_ev, u,v,w, derivs, u!=end_[0], v!=end_[1], w!=end_[2]);
		for(uint j=0; j<pts.size(); j++)
			pts[j] += basis_ev[j]*cp;
	}
}

void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w, Go::BasisDerivs2 & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	std::vector<double> values;
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.prepareDerivs(param_u, param_v, param_w, 0, 0, 0, nPts);

	int i=0;
	
	for(it=itStart; it!=itStop; ++it, ++i) {
		(*it)->evaluate(values, param_u, param_v, param_w, 2, param_u!=end_[0], param_v!=end_[1], param_w!=end_[2]);
	
		result.basisValues[i]    = values[0];
		result.basisDerivs_u[i]  = values[1];
		result.basisDerivs_v[i]  = values[2];
		result.basisDerivs_w[i]  = values[3];
		result.basisDerivs_uu[i] = values[4];
		result.basisDerivs_uv[i] = values[5];
		result.basisDerivs_uw[i] = values[6];
		result.basisDerivs_vv[i] = values[7];
		result.basisDerivs_vw[i] = values[8];
		result.basisDerivs_ww[i] = values[9];
	}
}

void LRSplineVolume::computeBasis (double param_u, double param_v, double param_w, Go::BasisDerivs & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	std::vector<double> values;
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.prepareDerivs(param_u, param_v, param_w, 0, 0, 0, nPts);
	
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i) {
		(*it)->evaluate(values, param_u, param_v, param_w, 1, param_u!=end_[0], param_v!=end_[1], param_w!=end_[2]);
		
		result.basisValues[i]   = values[0];
		result.basisDerivs_u[i] = values[1];
		result.basisDerivs_v[i] = values[2];
		result.basisDerivs_w[i] = values[3];
	}
}


void LRSplineVolume::computeBasis(double param_u, double param_v, double param_w, Go::BasisPts & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.preparePts(param_u, param_v, param_w, 0, 0, 0, nPts);
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i)
		result.basisValues[i] = (*it)->evaluate(param_u, param_v, param_w, param_u!=end_[0], param_v!=end_[1], param_w!=end_[2]);
}

void LRSplineVolume::computeBasis (double param_u,
                                   double param_v,
                                   double param_w,
                                   std::vector<std::vector<double> >& result,
                                   int derivs,
                                   int iEl ) const
{
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	result.clear();
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();
	result.resize(nPts);

	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i)
	    (*it)->evaluate(result[i], param_u, param_v, param_w, derivs, param_u!=end_[0], param_v!=end_[1], param_w!=end_[2]);
}

int LRSplineVolume::getElementContaining(double u, double v, double w) const {
	int ans = 0;
	for(Element *el : element_) {
		if(el->getParmin(0) <= u && el->getParmin(1) <= v && el->getParmin(2) <= w) 
			if((u < el->getParmax(0) || (u == end_[0] && u <= el->getParmax(0))) && 
			   (v < el->getParmax(1) || (v == end_[1] && v <= el->getParmax(1))) && 
			   (w < el->getParmax(2) || (w == end_[2] && w <= el->getParmax(2))))
				return ans;
		ans++;
	}
	return -1;
}

void LRSplineVolume::getMinspanRects(int iEl, std::vector<MeshRectangle*>& lines) {
	Element *e = element_[iEl];
	std::vector<Basisfunction*>::iterator it;
	double umin        = e->getParmin(0);
	double umax        = e->getParmax(0);
	double vmin        = e->getParmin(1);
	double vmax        = e->getParmax(1);
	double wmin        = e->getParmin(2);
	double wmax        = e->getParmax(2);
	double min_du      = DBL_MAX;
	double min_dv      = DBL_MAX;
	double min_dw      = DBL_MAX;
	int    best_startI = order_[0]+2;
	int    best_stopI  = order_[0]+2;
	int    best_startJ = order_[1]+2;
	int    best_stopJ  = order_[1]+2;
	int    best_startK = order_[2]+2;
	int    best_stopK  = order_[2]+2;
	double du          = umax - umin;
	double dv          = vmax - vmin;
	double dw          = wmax - wmin;
	double midU        = (umax + umin)/2.0;
	double midV        = (vmax + vmin)/2.0;
	double midW        = (wmax + wmin)/2.0;
	bool   drop_u_rect = du >= maxAspectRatio_*dv || du >= maxAspectRatio_*dw;
	bool   drop_v_rect = dv >= maxAspectRatio_*du || dv >= maxAspectRatio_*dw;
	bool   drop_w_rect = dw >= maxAspectRatio_*du || dw >= maxAspectRatio_*dv;
	// loop over all supported B-splines and choose the minimum one
	for(Basisfunction *b : e->support() ) {
		double lowu  = b->getParmin(0);
		double highu = b->getParmax(0);
		double lowv  = b->getParmin(1);
		double highv = b->getParmax(1);
		double loww  = b->getParmin(2);
		double highw = b->getParmax(2);
		du = highu - lowu;
		dv = highv - lowv;
		dw = highv - lowv;
		int startI=0;
		int stopI=0;
		int startJ=0;
		int stopJ=0;
		int startK=0;
		int stopK=0;
		while((*b)[0][startI] <= e->getParmin(0))
			startI++;
		while((*b)[0][stopI]  <  e->getParmax(0))
			stopI++;
		while((*b)[1][startJ] <= e->getParmin(1))
			startJ++;
		while((*b)[1][stopJ]  <  e->getParmax(1))
			stopJ++;
		while((*b)[2][startK] <= e->getParmin(2))
			startK++;
		while((*b)[2][stopK]  <  e->getParmax(2))
			stopK++;

		// min_du is defined as the minimum TOTAL knot span (of an entire basis function)
		bool fixU = false;
		int delta_startI = abs(startI - (order_[0]+1)/2);
		int delta_stopI  = abs(stopI  - (order_[0]+1)/2);
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
		int delta_startJ = abs(startJ - (order_[1]+1)/2);
		int delta_stopJ  = abs(stopJ  - (order_[1]+1)/2);
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

		bool fixW = false;
		int delta_startK = abs(startK - (order_[2]+1)/2);
		int delta_stopK  = abs(stopK  - (order_[2]+1)/2);
		if(  dw <  min_dw )
			fixW = true;
		if( dw == min_dw && delta_startK <= best_startK && delta_stopK  <= best_stopK )
			fixW = true;
		if(fixW) {
			wmin = loww;
			wmax = highw;
			min_dw = wmax-wmin;
			best_startK = delta_startK;
			best_stopK  = delta_stopK;
		}
	}
	if(! drop_u_rect )
		lines.push_back(new MeshRectangle(midU, vmin, wmin,  midU, vmax, wmax));
	if(! drop_v_rect )
		lines.push_back(new MeshRectangle(umin, midV, wmin,  umax, midV, wmax));
	if(! drop_w_rect )
		lines.push_back(new MeshRectangle(umin, vmin, midW,  umax, vmax, midW));
}

void LRSplineVolume::getFullspanRects(int iEl, std::vector<MeshRectangle*>& lines) {
	std::vector<Basisfunction*>::iterator it;
	Element *e = element_[iEl];
	double umin        = e->getParmin(0);
	double umax        = e->getParmax(0);
	double vmin        = e->getParmin(1);
	double vmax        = e->getParmax(1);
	double wmin        = e->getParmin(2);
	double wmax        = e->getParmax(2);
	double du          = umax - umin;
	double dv          = vmax - vmin;
	double dw          = wmax - wmin;
	double midU        = (umax + umin)/2.0;
	double midV        = (vmax + vmin)/2.0;
	double midW        = (wmax + wmin)/2.0;
	bool   drop_u_rect = du >= maxAspectRatio_*dv || du >= maxAspectRatio_*dw;
	bool   drop_v_rect = dv >= maxAspectRatio_*du || dv >= maxAspectRatio_*dw;
	bool   drop_w_rect = dw >= maxAspectRatio_*du || dw >= maxAspectRatio_*dv;
	// loop over all supported B-splines and make sure that everyone is covered by meshrect
	for(Basisfunction *b : e->support() ) {
		umin = (umin > b->getParmin(0)) ? b->getParmin(0) : umin;
		umax = (umax < b->getParmax(0)) ? b->getParmax(0) : umax;
		vmin = (vmin > b->getParmin(1)) ? b->getParmin(1) : vmin;
		vmax = (vmax < b->getParmax(1)) ? b->getParmax(1) : vmax;
		wmin = (wmin > b->getParmin(2)) ? b->getParmin(2) : wmin;
		wmax = (wmax < b->getParmax(2)) ? b->getParmax(2) : wmax;
	}

	if(! drop_u_rect )
		lines.push_back(new MeshRectangle(midU, vmin, wmin,  midU, vmax, wmax));
	if(! drop_v_rect )
		lines.push_back(new MeshRectangle(umin, midV, wmin,  umax, midV, wmax));
	if(! drop_w_rect )
		lines.push_back(new MeshRectangle(umin, vmin, midW,  umax, vmax, midW));
}

void LRSplineVolume::getStructMeshRects(Basisfunction *b, std::vector<MeshRectangle*>& rects) {
	double umin = b->getParmin(0);
	double umax = b->getParmax(0);
	double vmin = b->getParmin(1);
	double vmax = b->getParmax(1);
	double wmin = b->getParmin(2);
	double wmax = b->getParmax(2);

	// find the largest knotspan in this function
	double max[3] = {0,0,0};
	for(int d=0; d<3; d++) {
		for(int j=0; j<order_[d]; j++) {
			double du = (*b)[d][j+1]-(*b)[d][j];
			bool isZeroSpan =  fabs(du) < DOUBLE_TOL ;
			max[d] = (isZeroSpan || max[d]>du) ? max[d] : du;
		}
	}

	// to keep as "square" basis function as possible, only insert
	// into the largest knot spans
	for(int j=0; j<order_[0]; j++) {
		double du = (*b)[0][j+1]-(*b)[0][j];
		if( fabs(du-max[0]) < DOUBLE_TOL )
			rects.push_back(new MeshRectangle(((*b)[0][j] + (*b)[0][j+1])/2.0, vmin, wmin,
			                                  ((*b)[0][j] + (*b)[0][j+1])/2.0, vmax, wmax));
	}
	for(int j=0; j<order_[1]; j++) {
		double dv = (*b)[1][j+1]-(*b)[1][j];
		if( fabs(dv-max[1]) < DOUBLE_TOL )
			rects.push_back(new MeshRectangle(umin, ((*b)[1][j] + (*b)[1][j+1])/2.0, wmin,
			                                  umax, ((*b)[1][j] + (*b)[1][j+1])/2.0, wmax));
	}
	for(int j=0; j<order_[2]; j++) {
		double dw = (*b)[2][j+1]-(*b)[2][j];
		if( fabs(dw-max[2]) < DOUBLE_TOL )
			rects.push_back(new MeshRectangle(umin, vmin, ((*b)[2][j] + (*b)[2][j+1])/2.0,
			                                  umax, vmax, ((*b)[2][j] + (*b)[2][j+1])/2.0));
	}
}

void LRSplineVolume::refineBasisFunction(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineBasisFunction(tmp);
}

void LRSplineVolume::refineBasisFunction(const std::vector<int> &indices) {
	std::vector<int> sortedInd(indices);
	std::sort(sortedInd.begin(), sortedInd.end());
	std::vector<MeshRectangle*> newRects;

	/* first retrieve all meshrects needed */
	int ib = 0;
	HashSet_iterator<Basisfunction*> it = basis_.begin();
	for(uint i=0; i<sortedInd.size(); i++) {
		while(ib < sortedInd[i]) {
			++ib;
			++it;
		}
		getStructMeshRects(*it,newRects);
	}

	/* Do the actual refinement */
	for(MeshRectangle *m : newRects) 
		insert_line(m);

	/* do a posteriori fixes to ensure a proper mesh */
	// aPosterioriFixes();

}

void LRSplineVolume::refineElement(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineElement(tmp);
}

void LRSplineVolume::refineElement(const std::vector<int> &indices) {
	std::vector<MeshRectangle*> newRects;

	/* first retrieve all meshrects needed */
	for(uint i=0; i<indices.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanRects(indices[i],newRects);
		else
			getFullspanRects(indices[i],newRects);
	}

	/* Do the actual refinement */
	for(uint i=0; i<newRects.size(); i++)
		insert_line(newRects[i]);

	/* do a posteriori fixes to ensure a proper mesh */
	// aPosterioriFixes();

}

void LRSplineVolume::refineByDimensionIncrease(const std::vector<double> &errPerElement, double beta) {
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
		}
		i++;
	} else {
		for(uint i=0; i<element_.size(); i++) 
			errors.push_back(IndexDouble(errPerElement[i], i));
	}

	/* sort errors */
	std::sort(errors.begin(), errors.end(), std::greater<IndexDouble>());

	/* first retrieve all possible meshrects needed */
	std::vector<std::vector<MeshRectangle*> > newRects(errors.size(), std::vector<MeshRectangle*>(0));
	for(uint i=0; i<errors.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanRects(errors[i].second, newRects[i]);
		else if(refStrat_ == LR_SAFE) 
			getFullspanRects(errors[i].second, newRects[i]);
		else if(refStrat_ == LR_ISOTROPIC_FUNC)  {
			Basisfunction *b = getBasisfunction(errors[i].second);
			getStructMeshRects(b, newRects[i]);
		}
		// note that this is an excessive loop as it computes the meshrects for ALL elements,
		// but we're only going to use a small part of this.
	}

	/* Do the actual refinement */
	int target_n_functions = ceil(basis_.size()*(1+beta));
	int i=0;
	while( basis_.size() < target_n_functions ) {
		for(uint j=0; j<newRects[i].size(); j++) {
			insert_line(newRects[i][j]);
		}
		i++;
	}

	/* do a posteriori fixes to ensure a proper mesh */
	// aPosterioriFixes();
}

#if 0
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
					if(fabs(left[j] - (vmin+vmax)/2) < best) {
						best = fabs(left[j] - (vmin+vmax)/2);
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
					if(fabs(right[j] - (vmin+vmax)/2) < best) {
						best = fabs(right[j] - (vmin+vmax)/2);
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
					if(fabs(top[j] - (umin+umax)/2) < best) {
						best = fabs(top[j] - (umin+umax)/2);
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
					if(fabs(bottom[j] - (umin+umax)/2) < best) {
						best = fabs(bottom[j] - (umin+umax)/2);
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
					getMinspanRects(i, splitLines);
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
#endif

MeshRectangle* LRSplineVolume::insert_line(MeshRectangle *newRect) {
	if(newRect->start_[0] < start_[0] ||
	   newRect->start_[1] < start_[1] ||
	   newRect->start_[2] < start_[2] ||
	   newRect->stop_[0]  > end_[0]  ||
	   newRect->stop_[1]  > end_[1]  ||
	   newRect->stop_[2]  > end_[2] ) {
		std::cerr << "Error: inserting meshrctangle " << *newRect << " outside parametric domain";
		std::cerr << "(" << start_[0] << ", " << start_[1] << ", " << start_[2] << ") x ";
		std::cerr << "(" <<   end_[0] << ", " <<   end_[1] << ", " <<   end_[2] << ") ";
		return NULL;
	}
	
	std::vector<MeshRectangle*> newGuys;
	newGuys.push_back(newRect);

	{ // check if the line is an extension or a merging of existing lines
#ifdef TIME_LRSPLINE
	PROFILE("meshrectangle verification");
#endif
	for(uint i=0; i<meshrect_.size(); i++) {
		for(uint j=0; j<newGuys.size(); j++) {
			int status = meshrect_[i]->makeOverlappingRects(newGuys, j);
			if(status == 1) { // deleted j, i kept unchanged
				j--;
				continue;
			} else if(status == 2) { // j kept unchanged, delete i
				delete meshrect_[i];
				meshrect_.erase(meshrect_.begin() + i);
				i--;
				break;
			} else if(status == 3) { // j kept unchanged, i added to newGuys
				meshrect_.erase(meshrect_.begin() + i);
				i--;
				break;
			} else if(status == 4) { // deleted j, i added to newGuys
				meshrect_.erase(meshrect_.begin() + i);
				i--;
				j--;
				break;
			}
		}
	}
	bool change = true;
	while(change) {
		change = false;
		for(uint i=0; i<newGuys.size() ; i++) {
			for(uint j=i+1; j<newGuys.size() && j>i; j++) {
				int status = newGuys[i]->makeOverlappingRects(newGuys, j);
				if(status == 1) { //deleted j, i kept unchanged
					j--;
				} else if(status == 2) { // j kept unchanged, delete i
					delete newGuys[i];
					newGuys.erase(newGuys.begin() + i);
					i--;
				} else if(status == 3) { // j kept unchanged, i added to newGuys
					newGuys.erase(newGuys.begin() + i);
					i--;
				} else if(status == 4) { // deleted j, i added to newGuys
					newGuys.erase(newGuys.begin() + i);
					i--;
				}
				if(status > 0) {
					change = true;
				}
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
		for(MeshRectangle *m : newGuys) {
			if(m->splits(b)) {
				int nKnots;
				if( ((nKnots=m->nKnotsIn(b)) != m->multiplicity_) ) {
					removeFunc.insert(b);
					split( m->constDirection(), b, m->constParameter(), m->multiplicity_-nKnots, newFuncStp1 ); 
				}
			}
		}
	}
	for(Basisfunction* b : removeFunc) {
		basis_.erase(b);
		delete b;
	}
	for(uint i=0; i<element_.size(); i++) {
		for(MeshRectangle *m : newGuys) {
			if(m->splits(element_[i]))
				element_.push_back(element_[i]->split(m->constDirection(), m->constParameter()) );
		}
	}
	} // end step 1 timer

	{ // STEP 2: test every NEW function against ALL old meshrects
#ifdef TIME_LRSPLINE
	PROFILE("STEP 2");
#endif
	for(MeshRectangle *m : newGuys) 
		meshrect_.push_back(m);
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
	return NULL;
}

void LRSplineVolume::split(int constDir, Basisfunction *b, double new_knot, int multiplicity, HashSet<Basisfunction*> &newFunctions) {
#ifdef TIME_LRSPLINE
	PROFILE("split()");
#endif

	// create the new functions b1 and b2
	Basisfunction *b1, *b2;
	std::vector<double> knot = b->getknots(constDir);
	int     p                = b->getOrder(constDir);
	int     insert_index     = 0;
	if(new_knot < knot[0] || knot[p] < new_knot)
		return;
	while(knot[insert_index] < new_knot)
		insert_index++;
	double alpha1 = (insert_index == p)  ? 1.0 : (new_knot-knot[0])/(knot[p-1]-knot[0]);
	double alpha2 = (insert_index == 1 ) ? 1.0 : (knot[p]-new_knot)/(knot[p]-knot[1]);
	double newKnot[p+2];
	std::copy(knot.begin(), knot.begin()+p+1, newKnot+1);
	newKnot[0] = new_knot;
	std::sort(newKnot, newKnot + p+2);
	if(constDir == 0) {
		b1 = new Basisfunction(newKnot  ,  (*b)[1].begin(),  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha1);
		b2 = new Basisfunction(newKnot+1,  (*b)[1].begin(),  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha2);
	} else if(constDir == 1) {
		b1 = new Basisfunction((*b)[0].begin(), newKnot   ,  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), newKnot+1 ,  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha2);
	} else { // insert in w
		b1 = new Basisfunction((*b)[0].begin(), (*b)[1].begin(),  newKnot   , b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), (*b)[1].begin(),  newKnot+1 , b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha2);
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
	controlPoints.resize(order_[0]*order_[1]*order_[2]*dim_, 0);
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
		while(true) {
			int p    = order_[0]-1;
			int newI = -1;
			double z;
			if(       knotU.size() < (uint) startU+order_[0]   || knotU[startU+  order_[0]-1] != min) {
				z    = min;
				newI = startU;
			} else if(knotU.size() < (uint) startU+2*order_[0] || knotU[startU+2*order_[0]-1] != max ) {
				z    = max;
				newI = startU + order_[0];
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
		while(true) {
			int p    = order_[1]-1;
			int newI = -1;
			double z;
			if(       knotV.size() < (uint) startV+order_[1]   || knotV[startV+  order_[1]-1] != min) {
				z    = min;
				newI = startV;
			} else if(knotV.size() < (uint) startV+2*order_[1] || knotV[startV+2*order_[1]-1] != max ) {
				z = max;
				newI = startV + order_[1];
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
		while(true) {
			int p    = order_[2]-1;
			int newI = -1;
			double z;
			if(       knotW.size() < (uint) startW+order_[2]   || knotW[startW+  order_[2]-1] != min) {
				z    = min;
				newI = startW;
			} else if(knotW.size() < (uint) startW+2*order_[2] || knotW[startW+2*order_[2]-1] != max ) {
				z = max;
				newI = startW + order_[2];
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
		for(int w=startW; w<startW+order_[2]; w++)
			for(int v=startV; v<startV+order_[1]; v++)
				for(int u=startU; u<startU+order_[0]; u++)
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
	int n1 = knots_u.size() - order_[0];
	int n2 = knots_v.size() - order_[1];
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
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_[0]+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_[1]+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_[0]; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_[1]; j++) {
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
					int p = order_[0]-1;
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
					int p = order_[1]-1;
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
	int n1 = knots_u.size() - order_[0];
	int n2 = knots_v.size() - order_[1];
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
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_[0]+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_[1]+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_[0]; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_[1]; j++) {
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
					int p = order_[0]-1;
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
					int p = order_[1]-1;
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
	int n1 = knots_u.size() - order_[0];
	int n2 = knots_v.size() - order_[1];
	int fullDim = n1*n2;
	bool fullVerbose   = fullDim < 30  && nmb_bas < 50;
	bool sparseVerbose = fullDim < 250 && nmb_bas < 100;

	std::vector<std::vector<double > > C;  // rational projection matrix 

	for (int i = 0; i < nmb_bas; ++i) {
		int startU, startV;
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_[0]+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_[1]+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_[0]; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_[1]; j++) {
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
					int p = order_[0]-1;
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
					int p = order_[1]-1;
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
			if(fabs(C[j][i+zeroColumns]) > maxPivot) {
				maxPivot = fabs(C[j][i+zeroColumns]);
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
		for(int j=0; j<order_[0]+1; j++)
			b->knot_u_[j] = floor(b->knot_u_[j]/scale + 0.5);
		for(int j=0; j<order_[1]+1; j++)
			b->knot_v_[j] = floor(b->knot_v_[j]/scale + 0.5);
	}
	
	// scale all LRSplineVolume values
	start_[0] = floor(start_[0]/scale + 0.5);
	start_[1] = floor(start_[1]/scale + 0.5);
	end_[0]   = floor(end_[0]  /scale + 0.5);
	end_[1]   = floor(end_[1]  /scale + 0.5);

	return scale;
}
#endif

void LRSplineVolume::getDiagonalElements(std::vector<int> &result) const  {
	result.clear();
	for(uint i=0; i<element_.size(); i++) 
		if(element_[i]->getParmin(0) == element_[i]->getParmin(1) && element_[i]->getParmin(0) == element_[i]->getParmin(2))
			result.push_back(i);
}

void LRSplineVolume::getDiagonalBasisfunctions(std::vector<Basisfunction*> &result) const  {
	result.clear();
	for(Basisfunction *b : basis_) {
		bool isDiag = true;
		for(int j=0; j<order_[0]+1; j++)
			if(b->getknots(0)[j] != b->getknots(1)[j] || b->getknots(0)[j] != b->getknots(2)[j])
				isDiag = false;
		if(isDiag)
			result.push_back(b);
	}
}


void LRSplineVolume::read(std::istream &is) {
	start_[0] =  DBL_MAX;
	end_[0]   = -DBL_MAX;
	start_[1] =  DBL_MAX;
	end_[1]   = -DBL_MAX;
	start_[2] =  DBL_MAX;
	end_[2]   = -DBL_MAX;

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
	is >> order_[0];        ws(is);
	is >> order_[1];        ws(is);
	is >> order_[2];        ws(is);
	is >> nBasis;          ws(is);
	is >> nMeshRectangles; ws(is);
	is >> nElements;       ws(is);
	is >> dim_;            ws(is);
	is >> rational_;       ws(is);
	
	meshrect_.resize(nMeshRectangles);
	element_.resize(nElements);
	basisVector.resize(nBasis);
	int allOrder[] = {order_[0], order_[1], order_[2]};

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	// read all basisfunctions
	for(int i=0; i<nBasis; i++) {
		Basisfunction *b = new Basisfunction(dim_, 3, allOrder);
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
		start_[0] = (element_[i]->getParmin(0) < start_[0]) ? element_[i]->getParmin(0) : start_[0];
		end_[0]   = (element_[i]->getParmax(0) > end_[0]  ) ? element_[i]->getParmax(0) : end_[0]  ;
		start_[1] = (element_[i]->getParmin(1) < start_[1]) ? element_[i]->getParmin(1) : start_[1];
		end_[1]   = (element_[i]->getParmax(1) > end_[1]  ) ? element_[i]->getParmax(1) : end_[1]  ;
		start_[2] = (element_[i]->getParmin(2) < start_[2]) ? element_[i]->getParmin(2) : start_[2];
		end_[2]   = (element_[i]->getParmax(2) > end_[2]  ) ? element_[i]->getParmax(2) : end_[2]  ;
	}
}

void LRSplineVolume::write(std::ostream &os) const {
	generateIDs();
	os << std::setprecision(16);
	os << "# LRSPLINE VOLUME\n";
	os << "#\tp1\tp2\tp3\tNbasis\tNline\tNel\tdim\trat\n\t";
	os << order_[0] << "\t";
	os << order_[1] << "\t";
	os << order_[2] << "\t";
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

#undef DOUBLE_TOL

} // end namespace LR

