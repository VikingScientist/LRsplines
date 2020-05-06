#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/MeshRectangle.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"

#include <algorithm>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
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

	initMeta();
}

#ifdef HAS_GOTOOLS
LRSplineVolume::LRSplineVolume(Go::SplineVolume *vol) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif

	initMeta();
	initCore(vol->numCoefs(0),      vol->numCoefs(1),      vol->numCoefs(2),
	         vol->order(0),         vol->order(1),         vol->order(2),
	         vol->basis(0).begin(), vol->basis(1).begin(), vol->basis(2).begin(),
			 vol->ctrl_begin(),     vol->dimension(),      vol->rational() );
}
#endif

/*
LRSplineVolume::LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, double *knot_u, double *knot_v, double *knot_w, double *coef, int dim, bool rational) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();
	initCore(n1, n2, n3, order_u, order_v, order_w, knot_u, knot_v, knot_w, coef, dim, rational);
}
*/

/************************************************************************************************************************//**
 * \brief Constructor. Creates a tensor product LRSplineVolume of the 3D unit cube
 * \param n1 The number of basisfunctions in the first parametric direction
 * \param n2 The number of basisfunctions in the second parametric direction
 * \param n3 The number of basisfunctions in the third parametric direction
 * \param order_u The order (polynomial degree + 1) in the first parametric direction
 * \param order_v The order (polynomial degree + 1) in the second parametric direction
 * \param order_w The order (polynomial degree + 1) in the third parametric direction
 * \param knot_u The first knot vector consisting of n1+order_u elements
 * \param knot_v The second knot vector consisting of n2+order_v elements
 * \param knot_w The third knot vector consisting of n2+order_v elements
 ***************************************************************************************************************************/
LRSplineVolume::LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, double *knot_u, double *knot_v, double *knot_w) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();

	std::vector<double> grev_u = LRSpline::getGrevillePoints(order_u, knot_u, knot_u + n1 + order_u);
	std::vector<double> grev_v = LRSpline::getGrevillePoints(order_v, knot_v, knot_v + n2 + order_v);
	std::vector<double> grev_w = LRSpline::getGrevillePoints(order_w, knot_w, knot_w + n3 + order_w);
	std::vector<double> coef(n1*n2*n3*3);
	int l=0;
	for(int k=0; k<n3; k++)
		for(int j=0; j<n2; j++)
			for(int i=0; i<n1; i++) {
				coef[l++] = grev_u[i];
				coef[l++] = grev_v[j];
				coef[l++] = grev_w[k];
	}

	initCore(n1, n2, n3, order_u, order_v, order_w, knot_u, knot_v, knot_w, coef.begin(), 3, false);
}

/************************************************************************************************************************//**
 * \brief Constructor. Creates a uniform tensor product LRSplineVolume of the 2D unit square
 * \param n1 The number of basisfunctions in the first parametric direction
 * \param n2 The number of basisfunctions in the second parametric direction
 * \param n3 The number of basisfunctions in the third parametric direction
 * \param order_u The order (polynomial degree + 1) in the first parametric direction
 * \param order_v The order (polynomial degree + 1) in the second parametric direction
 * \param order_w The order (polynomial degree + 1) in the third parametric direction
 ***************************************************************************************************************************/
LRSplineVolume::LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();
	// sanity check input
	if(n1 < order_u ||
	   n2 < order_v ||
	   n3 < order_w) {
		std::cerr << "Error: n<p in LRSplineVolume constructor\n";
		// really ought to throw exception here, but havent the framework
		// for this up and running yet. Make it a zombie surface
		double knot[4] = {0,0,1,1};
		double cp[4] = {0,0,0,0};
		initCore(2,2,2,2,2,2,knot,knot,knot,cp,1); // init dummy-state object
		return;
	}

	std::vector<double> knot_u = LRSpline::getUniformKnotVector(n1, order_u);
	std::vector<double> knot_v = LRSpline::getUniformKnotVector(n2, order_v);
	std::vector<double> knot_w = LRSpline::getUniformKnotVector(n3, order_w);
	std::vector<double> grev_u = LRSpline::getGrevillePoints(order_u, knot_u.begin(), knot_u.end());
	std::vector<double> grev_v = LRSpline::getGrevillePoints(order_v, knot_v.begin(), knot_v.end());
	std::vector<double> grev_w = LRSpline::getGrevillePoints(order_w, knot_w.begin(), knot_w.end());
	std::vector<double> coef(n1*n2*n3*3);
	int l=0;
	for(int k=0; k<n3; k++)
		for(int j=0; j<n2; j++)
			for(int i=0; i<n1; i++) {
				coef[l++] = grev_u[i];
				coef[l++] = grev_v[j];
				coef[l++] = grev_w[k];
	}

	// generate the uniform knot vector
	initCore(n1,n2,n3, order_u, order_v, order_w, knot_u.begin(), knot_v.begin(), knot_w.begin(), coef.begin(), 3, false);
}

LRSplineVolume::~LRSplineVolume() {
	for(Basisfunction* b : basis_)
		delete b;
	for(uint i=0; i<meshrect_.size(); i++)
		delete meshrect_[i];
	for(uint i=0; i<element_.size(); i++)
		delete element_[i];
}

void LRSplineVolume::initMeta() {
	maxTjoints_           = -1;
	doCloseGaps_          = true;
	maxAspectRatio_       = 2.0;
	doAspectRatioFix_     = false;
	refStrat_             = LR_FULLSPAN;
	refKnotlineMult_      = 1;
	symmetry_             = 1;
	builtElementCache_    = false;
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

#ifdef HAS_GOTOOLS
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
	if(iEl == -1)
		return;

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
	if(iEl == -1)
		return;

	for(Basisfunction* b : element_[iEl]->support() ) {
		b->getControlPoint(cp);
		b->evaluate(basis_ev, u,v,w, derivs, u!=end_[0], v!=end_[1], w!=end_[2]);
		for(uint j=0; j<pts.size(); j++)
			pts[j] += basis_ev[j]*cp;
	}
}
#endif

/************************************************************************************************************************//**
 * \brief Evaluate the volume at a point (u,v,w)
 * \param[out] pt The result, i.e. the parametric volume mapped to physical space
 * \param u The u-coordinate on which to evaluate the volume
 * \param v The v-coordinate on which to evaluate the volume
 * \param w The w-coordinate on which to evaluate the volume
 * \param iEl The element index which this point is contained in. If used will speed up computational efficiency
 * \param u_from_right True if first coordinate should be evaluated in the limit from the right
 * \param v_from_right True if second coordinate should be evaluated in the limit from the right
 * \param w_from_right True if third coordinate should be evaluated in the limit from the right
 ***************************************************************************************************************************/
void LRSplineVolume::point(std::vector<double> &pt, double u, double v, double w, int iEl, bool u_from_right, bool v_from_right, bool w_from_right) const {
	std::vector<std::vector<double> > res;
	point(res, u, v, w, 0, u_from_right, v_from_right, w_from_right, iEl);
	pt = res[0];
}

/************************************************************************************************************************//**
 * \brief Evaluate the volume at a point (u,v)
 * \param[out] pt The result, i.e. the parametric volume mapped to physical space
 * \param u The u-coordinate on which to evaluate the volume
 * \param v The v-coordinate on which to evaluate the volume
 * \param w The w-coordinate on which to evaluate the volume
 * \param iEl The element index which this point is contained in. If used will speed up computational efficiency
 ***************************************************************************************************************************/
void LRSplineVolume::point(std::vector<double> &pt, double u, double v, double w, int iEl) const {
	std::vector<std::vector<double> > res;
	point(res, u, v, w, 0, iEl);
	pt = res[0];
}

/************************************************************************************************************************//**
 * \brief Evaluate the volume and its derivatives at a point (u,v)
 * \param[out] pts The result, i.e. the parametric volume as well as all parametric derivatives
 * \param u The u-coordinate on which to evaluate the volume
 * \param v The v-coordinate on which to evaluate the volume
 * \param w The w-coordinate on which to evaluate the volume
 * \param derivs The number of derivatives requested
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details Pts is a vector where the first element is the values themselves, the next two is the first derivatives du, dv and dw
 *          The next six is the second derivatives in order d2u, dudv, dudw, d2v, dvdw and d2w. The next ten are the third derivatives
 ***************************************************************************************************************************/
void LRSplineVolume::point(std::vector<std::vector<double> > &pts, double u, double v, double w, int derivs, int iEl) const {
	point(pts, u, v, w, derivs, u!=end_[0], v!=end_[1], w!=end_[2], iEl);
}

/************************************************************************************************************************//**
 * \brief Evaluate the volume and its derivatives at a point (u,v)
 * \param[out] pts The result, i.e. the parametric volume as well as all parametric derivatives
 * \param u The u-coordinate on which to evaluate the volume
 * \param v The v-coordinate on which to evaluate the volume
 * \param w The w-coordinate on which to evaluate the volume
 * \param derivs The number of derivatives requested
 * \param u_from_right True if first coordinate should be evaluated in the limit from the right
 * \param v_from_right True if second coordinate should be evaluated in the limit from the right
 * \param w_from_right True if third coordinate should be evaluated in the limit from the right
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details Pts is a vector where the first element is the values themselves, the next three is the first derivatives du, dv and dw
 *          The next six is the second derivatives in order d2u, dudv, dudw, d2v, dvdw and d2w. The next ten are the third derivatives
 ***************************************************************************************************************************/
void LRSplineVolume::point(std::vector<std::vector<double> > &pts, double u, double v, double w, int derivs, bool u_from_right, bool v_from_right, bool w_from_right, int iEl) const {
#ifdef TIME_LRSPLINE
	PROFILE("Point()");
#endif
	std::vector<double> basis_ev;

	// clear and resize output array (optimization may consider this an outside task)
	pts.clear();
	pts.resize((derivs+1)*(derivs+2)*(derivs+3)/6);
	for(uint i=0; i<pts.size(); i++)
		pts[i].resize(dim_, 0);

	if(iEl == -1)
		iEl = getElementContaining(u,v,w);
	if(iEl == -1)
		return;
	for(Basisfunction* b : element_[iEl]->support() ) {
		b->evaluate(basis_ev, u,v,w, derivs, u_from_right, v_from_right, w_from_right);
		for(uint i=0; i<pts.size(); i++)
			for(int j=0; j<dim_; j++)
				pts[i][j] += basis_ev[i]*b->cp(j);
	}
}

#ifdef HAS_GOTOOLS
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
#endif

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

/************************************************************************************************************************//**
 * \brief Computes a cached lookup table for quick determination of element distribution. Allows getElementContaining()
 *        to be ran in O(log(n)) time.
 * \details This is done by creating a full tensor mesh (of elements) and storing the LR-elements from which each
 *          sub-element came from.
 ***************************************************************************************************************************/
void LRSplineVolume::createElementCache() const {
	generateIDs();
	// find the set of all unique knots in each direction (this is our global mesh)
	glob_knot_u_.clear();
	glob_knot_v_.clear();
	glob_knot_w_.clear();
	this->getGlobalUniqueKnotVector(glob_knot_u_, glob_knot_v_, glob_knot_w_);

	// create a tensor-mesh of elements given by an nxm matrix
	elementCache_ = std::vector<std::vector<std::vector<int> > >(glob_knot_u_.size(),
                              std::vector<std::vector<int> >(  glob_knot_v_.size(),
                                          std::vector<int>(    glob_knot_w_.size(), -1)));

	// store parent element (i.e. the actual LR element) for each sub-cell
	for(Element *e : getAllElements()) {
		// binary-search to look up indices
		int i0 = std::lower_bound(glob_knot_u_.begin(), glob_knot_u_.end(), e->umin()) - glob_knot_u_.begin();
		int i1 = std::lower_bound(glob_knot_u_.begin(), glob_knot_u_.end(), e->umax()) - glob_knot_u_.begin();
		int j0 = std::lower_bound(glob_knot_v_.begin(), glob_knot_v_.end(), e->vmin()) - glob_knot_v_.begin();
		int j1 = std::lower_bound(glob_knot_v_.begin(), glob_knot_v_.end(), e->vmax()) - glob_knot_v_.begin();
		int k0 = std::lower_bound(glob_knot_w_.begin(), glob_knot_w_.end(), e->wmin()) - glob_knot_w_.begin();
		int k1 = std::lower_bound(glob_knot_w_.begin(), glob_knot_w_.end(), e->wmax()) - glob_knot_w_.begin();
		for(int i=i0; i<i1; i++)
			for(int j=j0; j<j1; j++)
			  for(int k=k0; k<k1; k++)
				  elementCache_[i][j][k] = e->getId();
	}
	builtElementCache_ = true;
}

/************************************************************************************************************************//**
 * \brief Get the element index of the element containing the parametric point (u,v,w)
 * \param u The u-coordinate
 * \param v The v-coordinate
 * \param w The w-coordinate
 * \return The index of the element which contains (u,v,w)
 ***************************************************************************************************************************/
int LRSplineVolume::getElementContaining(double u, double v, double w) const {
	// sanity check input
	if(u < startparam(0) || u > endparam(0) || v < startparam(1) || v > endparam(1) || w < startparam(2) || w > endparam(2))
		return -1;
	// build cache if not already present
	if(builtElementCache_ == false)
		createElementCache();

	// binary search for the right element
	size_t i = std::upper_bound(glob_knot_u_.begin(), glob_knot_u_.end(), u) - glob_knot_u_.begin() - 1;
	size_t j = std::upper_bound(glob_knot_v_.begin(), glob_knot_v_.end(), v) - glob_knot_v_.begin() - 1;
	size_t k = std::upper_bound(glob_knot_w_.begin(), glob_knot_w_.end(), w) - glob_knot_w_.begin() - 1;

	// special case endpoints since element definitions are defined to be [umin,umax) for all except last element: [umin,umax]
	if(i == glob_knot_u_.size()-1) i--;
	if(j == glob_knot_v_.size()-1) j--;
	if(k == glob_knot_w_.size()-1) k--;
	return elementCache_[i][j][k];
}

/************************************************************************************************************************//**
 * \brief Get the all neighbours for a given Element
 * \param iEl  The element index
 * \param edge The boundary edge
 * \return A set of indexes of all neighbouring elements on a given edge
 ***************************************************************************************************************************/
std::set<int> LRSplineVolume::getElementNeighbours(int iEl, parameterEdge edge) const {
	// build cache if not already present
	if(builtElementCache_ == false)
		createElementCache();

	Element *e = element_[iEl];
	int i0 = std::lower_bound(glob_knot_u_.begin(), glob_knot_u_.end(), e->umin()) - glob_knot_u_.begin();
	int i1 = std::lower_bound(glob_knot_u_.begin(), glob_knot_u_.end(), e->umax()) - glob_knot_u_.begin();
	int j0 = std::lower_bound(glob_knot_v_.begin(), glob_knot_v_.end(), e->vmin()) - glob_knot_v_.begin();
	int j1 = std::lower_bound(glob_knot_v_.begin(), glob_knot_v_.end(), e->vmax()) - glob_knot_v_.begin();
	int k0 = std::lower_bound(glob_knot_w_.begin(), glob_knot_w_.end(), e->wmin()) - glob_knot_w_.begin();
	int k1 = std::lower_bound(glob_knot_w_.begin(), glob_knot_w_.end(), e->wmax()) - glob_knot_w_.begin();
	std::set<int> result;
	if (edge & (TOP|BOTTOM)) {
		size_t k = (edge == TOP) ? k1:k0-1;
		if (k >= 0 && k < glob_knot_w_.size()-1)
			for(int i=i0; i<i1; i++) {
				for(int j=j0; j<j1; j++) {
					result.insert(elementCache_[i][j][k]);
				}
			}
	} else if (edge & (NORTH|SOUTH)) {
		size_t j = (edge == NORTH) ? j1:j0-1;
		if (j >= 0 && j < glob_knot_v_.size()-1)
			for(int i=i0; i<i1; i++) {
				for(int k=k0; k<k1; k++) {
					result.insert(elementCache_[i][j][k]);
				}
			}
	} else if (edge & (WEST|EAST)) {
		size_t i = (edge == EAST) ? i1:i0-1;
		if (i >= 0 && i < glob_knot_u_.size()-1)
			for(int j=j0; j<j1; j++) {
				for(int k=k0; k<k1; k++) {
					result.insert(elementCache_[i][j][k]);
				}
			}
	}
	return result;
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
		if( fabs(du-max[0]) < DOUBLE_TOL ) {
			MeshRectangle *m = new MeshRectangle(((*b)[0][j] + (*b)[0][j+1])/2.0, vmin, wmin,
			                                     ((*b)[0][j] + (*b)[0][j+1])/2.0, vmax, wmax);
			if(!MeshRectangle::addUniqueRect(rects, m))
				delete m;
		}
	}
	for(int j=0; j<order_[1]; j++) {
		double dv = (*b)[1][j+1]-(*b)[1][j];
		if( fabs(dv-max[1]) < DOUBLE_TOL ) {
			MeshRectangle *m = new MeshRectangle(umin, ((*b)[1][j] + (*b)[1][j+1])/2.0, wmin,
			                                     umax, ((*b)[1][j] + (*b)[1][j+1])/2.0, wmax);
			if(!MeshRectangle::addUniqueRect(rects, m))
				delete m;
		}
	}
	for(int j=0; j<order_[2]; j++) {
		double dw = (*b)[2][j+1]-(*b)[2][j];
		if( fabs(dw-max[2]) < DOUBLE_TOL ) {
			MeshRectangle *m = new MeshRectangle(umin, vmin, ((*b)[2][j] + (*b)[2][j+1])/2.0,
			                                     umax, vmax, ((*b)[2][j] + (*b)[2][j+1])/2.0);
			if(!MeshRectangle::addUniqueRect(rects, m))
				delete m;
		}
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
	if(refStrat_ == LR_STRUCTURED_MESH) { // error per-function
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

	/* first retrieve all possible meshrects needed */
	std::vector<std::vector<MeshRectangle*> > newRects(errors.size(), std::vector<MeshRectangle*>(0));
	for(uint i=0; i<errors.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanRects(errors[i].second, newRects[i]);
		else if(refStrat_ == LR_FULLSPAN)
			getFullspanRects(errors[i].second, newRects[i]);
		else if(refStrat_ == LR_STRUCTURED_MESH)  {
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

std::vector<Meshline*> LRSplineVolume::getEdgeKnots(parameterEdge edge, bool normalized) const {
	std::vector<Meshline*> results;

	// fetch global patch size
	double u0 = this->startparam(0);
	double u1 = this->endparam(0);
	double v0 = this->startparam(1);
	double v1 = this->endparam(1);
	double w0 = this->startparam(2);
	double w1 = this->endparam(2);

	for(auto rect : this->getAllMeshRectangles()) {
		if(edge == WEST && fabs(rect->start_[0] - u0) < DOUBLE_TOL  ||
		   edge == EAST && fabs(rect->stop_[0]  - u1) < DOUBLE_TOL) {
			if(rect->constDirection() == 1)
				results.push_back(new Meshline(false, rect->start_[1], rect->start_[2], rect->stop_[2], rect->multiplicity_));
			else if(rect->constDirection() == 2)
				results.push_back(new Meshline(true, rect->start_[2], rect->start_[1], rect->stop_[1], rect->multiplicity_));

		} else if(edge == SOUTH && fabs(rect->start_[1] - v0) < DOUBLE_TOL  ||
		          edge == NORTH && fabs(rect->stop_[1]  - v1) < DOUBLE_TOL) {
			if(rect->constDirection() == 0)
				results.push_back(new Meshline(false, rect->start_[0], rect->start_[2], rect->stop_[2], rect->multiplicity_));
			else if(rect->constDirection() == 2)
				results.push_back(new Meshline(true, rect->start_[2], rect->start_[0], rect->stop_[0], rect->multiplicity_));

		} else if(edge == BOTTOM && fabs(rect->start_[2] - w0) < DOUBLE_TOL  ||
		          edge == TOP    && fabs(rect->stop_[2]  - w1) < DOUBLE_TOL) {
			if(rect->constDirection() == 0)
				results.push_back(new Meshline(false, rect->start_[0], rect->start_[1], rect->stop_[1], rect->multiplicity_));
			else if(rect->constDirection() == 1)
				results.push_back(new Meshline(true, rect->start_[1], rect->start_[0], rect->stop_[0], rect->multiplicity_));
		}
	}

	if(normalized) {
		for(auto r : results) {
			switch(edge) {
				case WEST:
				case EAST:
					if(r->span_u_line_) {
						r->start_ = (r->start_ - v0 ) / (v1-v0);
						r->stop_  = (r->stop_  - v0 ) / (v1-v0);
					} else {
						r->start_ = (r->start_ - w0 ) / (w1-w0);
						r->stop_  = (r->stop_  - w0 ) / (w1-w0);
					}
					r->const_par_ = (r->const_par_ - u0 ) / (u1-u0);
					break;
				case NORTH:
				case SOUTH:
					if(r->span_u_line_) {
						r->start_ = (r->start_ - u0 ) / (u1-u0);
						r->stop_  = (r->stop_  - u0 ) / (u1-u0);
					} else {
						r->start_ = (r->start_ - w0 ) / (w1-w0);
						r->stop_  = (r->stop_  - w0 ) / (w1-w0);
					}
					r->const_par_ = (r->const_par_ - v0 ) / (v1-v0);
					break;
				case TOP:
				case BOTTOM:
					if(r->span_u_line_) {
						r->start_ = (r->start_ - u0 ) / (u1-u0);
						r->stop_  = (r->stop_  - u0 ) / (u1-u0);
					} else {
						r->start_ = (r->start_ - v0 ) / (v1-v0);
						r->stop_  = (r->stop_  - v0 ) / (v1-v0);
					}
					r->const_par_ = (r->const_par_ - w0 ) / (w1-w0);
					break;
				default:
					break;
			}
		}
	}
	return results;
}

/************************************************************************************************************************//**
 * \brief Creates matching discretization (mesh) across patch boundaries so they can be stiched together
 * \param edge      which of the 8 parameter edges (NORTH/SOUTH/EAST/WEST/TOP/BOTTOM) that should be refined
 * \param functions list of normalized functions (knot vectors in range (0,1)) that should appear on the boundary edge of this patch
 * \details This will refine this patche so it will have conforming mesh on the interface between this patch and another. It will
 *          not manipulate the control-points of the basis functions, and it is left to the user to make sure that these coincide
 *          if the mesh is to match in the physical space, and not just in the parametric space.
 ***************************************************************************************************************************/
void LRSplineVolume::matchParametricEdge(parameterEdge edge, const std::vector<Basisfunction*> &functions) {
	double u0 = this->startparam(0);
	double u1 = this->endparam(0);
	double v0 = this->startparam(1);
	double v1 = this->endparam(1);
	double w0 = this->startparam(2);
	double w1 = this->endparam(2);
	for(auto b : functions) {
		for(int d=0; d<3; d++) {
			int mult = 1;
			auto knots = b->getknots(d);
			for(uint i=0; i<knots.size(); i++) {
				if( i==knots.size()-1 || fabs(knots[i+1] - knots[i])>DOUBLE_TOL) {
					if(d==0) {
						this->insert_line(new MeshRectangle((u1+u0)*knots[i]+u0, (v1-v0)*b->getParmin(1)+v0, (w1-w0)*b->getParmin(2)+w0,
						                                    (u1+u0)*knots[i]+u0, (v1-v0)*b->getParmax(1)+v0, (w1-w0)*b->getParmax(2)+w0, mult));
					} else if(d==1) {
						this->insert_line(new MeshRectangle((u1+u0)*b->getParmin(0)+u0, (v1-v0)*knots[i]+v0, (w1-w0)*b->getParmin(2)+w0,
						                                    (u1+u0)*b->getParmax(0)+u0, (v1-v0)*knots[i]+v0, (w1-w0)*b->getParmax(2)+w0, mult));
					} else {
						this->insert_line(new MeshRectangle((u1+u0)*b->getParmin(0)+u0, (v1-v0)*b->getParmin(1)+v0, (w1-w0)*knots[i]+w0,
						                                    (u1+u0)*b->getParmax(0)+u0, (v1-v0)*b->getParmax(1)+v0, (w1-w0)*knots[i]+w0, mult));
					}
					mult = 1;
				} else {
					mult++;
				}
			}
		}
	}
	aPosterioriFixElements();
}

/************************************************************************************************************************//**
 * \brief Creates matching discretization (mesh) across patch boundaries so they can be stiched together in a C0-fashion
 * \param edge      which of the 4 parameter edges (NORTH/SOUTH/EAST/WEST) on the master patch (this patch) that should be matched
 * \param other     pointer to the other patch that is to be matched
 * \param otherEdge which of the 4 parameter edges on the slave patch that should be matched
 * \returns true if any new refinements were introduced
 * \details This will refine this patche so it will have conforming mesh on the interface between this patch and another. It will
 *          not manipulate the control-points of the basis functions, and it is left to the user to make sure that these coincide
 *          if the mesh is to match in the physical space, and not just in the parametric space.
 ***************************************************************************************************************************/
bool LRSplineVolume::matchParametricEdge(parameterEdge edge, LRSplineVolume *other, parameterEdge otherEdge, bool reverse_u, bool reverse_v, bool flip_uv) {
	int n1 = this->nBasisFunctions();
	int n2 = other->nBasisFunctions();
	std::vector<Basisfunction*> fun1, fun2, normalized1, normalized2;
	this->getEdgeFunctions(fun1, edge);
	other->getEdgeFunctions(fun2, otherEdge);

	// NOTE: normalizing is non-trivial when the matching boundaries are of different cross-parametric sizes. I.e. this example here:
	/*
	 *   +------+ +-----------------+  ^
	 *   |      | |                 |  |
	 *   |  #1  | |       #2        |  2
	 *   |      | |                 |  |
	 *   +------+ +-----------------+  v
	 *
	 *   <- 2 --> <------- 4 ------->
	 *
	 * Rescaling the functions of patch #2 in the u-direction to (0,1) may match a different scaling of patch #1 scaled in its u-direction
	 * Note that we DO need to match the u-components of the basifunctions even if this is a v-direction edge that is being matched!
	 */

	// make a pre-processed rescaled set of basisfunctions which will be passed to the other LRSplineSurface patch
	int thismax  = (     edge==EAST ||      edge==NORTH ||      edge==TOP);
	int othermax = (otherEdge==EAST || otherEdge==NORTH || otherEdge==TOP);

	// fetch local (u,v)-coordinates for this edge
	int u1,v1,w1,u2,v2,w2;
	switch(edge) {
		case WEST:
		case EAST:
			u1 = 1; // first running direction
			v1 = 2; // second running
			w1 = 0; // const direction
			break;
		case NORTH:
		case SOUTH:
			u1 = 0; // first running direction
			v1 = 2; // second running
			w1 = 1; // const direction
			break;
		case TOP:
		case BOTTOM:
			u1 = 0; // first running direction
			v1 = 1; // second running
			w1 = 2; // const direction
			break;
        default:
            return false; // do nothing for 'NONE' or corner-input (SOUTH_WEST)
	}
	switch(otherEdge) {
		case WEST:
		case EAST:
			u2 = 1;
			v2 = 2;
			w2 = 0;
			break;
		case NORTH:
		case SOUTH:
			u2 = 0;
			v2 = 2;
			w2 = 1;
			break;
		case TOP:
		case BOTTOM:
			u2 = 0;
			v2 = 1;
			w2 = 2;
			break;
        default:
            return false; // do nothing for 'NONE' or corner-input (SOUTH_WEST)
	}

	for(auto b : fun1) {
		normalized1.push_back(b->copy());
		normalized1.back()->normalize(0, this->startparam(0), this->endparam(0));
		normalized1.back()->normalize(1, this->startparam(1), this->endparam(1));
		normalized1.back()->normalize(2, this->startparam(2), this->endparam(2));
		if(flip_uv)
			normalized1.back()->flip(u1,v1);
		if(reverse_u)
			normalized1.back()->reverse(u1);
		if(reverse_v)
			normalized1.back()->reverse(v1);
		if(thismax != othermax)
			normalized1.back()->reverse(w1);
		normalized1.back()->flip(u1,u2);
		normalized1.back()->flip(v1,v2);
	}

	// make a pre-processed rescaled set of basisfunctions that will be passed from the other patch and onto this
	for(auto b : fun2) {
		normalized2.push_back(b->copy());
		normalized2.back()->normalize(0, other->startparam(0), other->endparam(0));
		normalized2.back()->normalize(1, other->startparam(1), other->endparam(1));
		normalized2.back()->normalize(2, other->startparam(2), other->endparam(2));
		if(flip_uv)
			normalized2.back()->flip(u2,v2);
		if(reverse_u)
			normalized2.back()->reverse(u2);
		if(reverse_v)
			normalized2.back()->reverse(v2);
		if(thismax != othermax)
			normalized2.back()->reverse(w2);
		normalized2.back()->flip(v1,v2);
		normalized2.back()->flip(u1,u2);
	}

	// refine patch #1
	this->matchParametricEdge(edge, normalized2);

	// refine patch #2
	other->matchParametricEdge(otherEdge, normalized1);

	// check if any change has happened:
	return (this->nBasisFunctions() != n1 || other->nBasisFunctions() != n2);

}


/************************************************************************************************************************//**
 * \brief Enforces all elements to be of equal size in all 3 directions
 * \return true if some elements were bisected
 * \details This is a rather specialized function and will only work in certain scenarios. The intended use is when 2:1-elements
 *          sneak into the mesh (for instance by using LRSplineSurface::matchParametricEdge()) and need to be further subdivided.
 *          They are then fixed by simpy inserting dividing this element in two along its longest direction. Unless this hits
 *          another meshline perfectly, this line WILL be too short for legal insertion and the method will not work.
 *
 *
 *          +-----------+
 *          |           |  <--- fix this element by splitting in two along the u-direction
 *          |           |
 *          +-----+-----+
 *          |     |     |
 *          |     |     |
 *          +-----+-----+
 *
 *                 \
 *                  \ won't work unless this particular line is here, such that the new line forms an elongation of the old one
 *
 *
 *
 ***************************************************************************************************************************/
bool LRSplineVolume::enforceIsotropic() {
	bool somethingFixed = true;
	bool result = false;
	while(somethingFixed) {
		somethingFixed = false;
		for(uint i=0; i<element_.size(); i++) {
			double umin = element_[i]->umin();
			double umax = element_[i]->umax();
			double vmin = element_[i]->vmin();
			double vmax = element_[i]->vmax();
			double wmin = element_[i]->wmin();
			double wmax = element_[i]->wmax();
			double du = umax-umin;
			double dv = vmax-vmin;
			double dw = wmax-wmin;
			double h  = std::min(std::min(du,dv),dw);
			MeshRectangle *m;
			if(du - h > DOUBLE_TOL) {
				m = new MeshRectangle(umin+du/2, vmin, wmin, umin+du/2, vmax, wmax, refKnotlineMult_);
				insert_line(m);
				somethingFixed = true;
			}
			if(dv - h > DOUBLE_TOL) {
				m = new MeshRectangle(umin, vmin+dv/2, wmin, umax, vmin+dv/2, wmax, refKnotlineMult_);
				insert_line(m);
				somethingFixed = true;
			}
			if(dw - h > DOUBLE_TOL) {
				m = new MeshRectangle(umin, vmin, wmin+dw/2, umax, vmax, wmin+dw/2, refKnotlineMult_);
				insert_line(m);
				somethingFixed = true;
			}
			if(somethingFixed) {
				result = true;
				break;
			}
		}
	}
	return result;
}

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
		exit(3239621);
	}
	if(newRect->multiplicity_ < 1) {
		std::cerr << "LRSplineVolume::insert_line() requested mesh rectangle with multiplicity " << newRect->multiplicity_ << ". Needs non-negative values\n";
		exit(5312174);
	}

	std::vector<MeshRectangle*> newGuys;
	newGuys.push_back(newRect);

	{ // check if the line is an extension or a merging of existing lines
#ifdef TIME_LRSPLINE
	PROFILE("meshrectangle verification");
#endif
	for(uint i=0; i<meshrect_.size(); i++) {
		for(uint j=0; j<newGuys.size(); j++) {
			int status = meshrect_[i]->makeOverlappingRects(newGuys, j, true);
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
			} else if(status == 5) { // deleted j, i duplicate in newGuys
				delete meshrect_[i];
				meshrect_.erase(meshrect_.begin() + i);
				i--;
				j--;
				break;
			} else if(status == 6) { // j kept unchanged, i duplicate in newGuys
				delete meshrect_[i];
				meshrect_.erase(meshrect_.begin() + i);
				i--;
				break;
			}
		}
	}
	bool change = true;
	while(change) {
		change = false;
		for(uint i=0; i<newGuys.size() && !change ; i++) {
			for(uint j=i+1; j<newGuys.size() && j>i; j++) {
				int status = newGuys[i]->makeOverlappingRects(newGuys, j, false);
				if(status == 1) { //deleted j, i kept unchanged
					;
				} else if(status == 2) { // j kept unchanged, deleted i
					delete newGuys[i];
					newGuys.erase(newGuys.begin() + i);
				} else if(status == 3) { // j kept unchanged, i added to newGuys
					newGuys.erase(newGuys.begin() + i);
				} else if(status == 4) { // deleted j, i added to newGuys
					newGuys.erase(newGuys.begin() + i);
				} else if(status == 5) { // deleted j, i duplicate in newGuys
					;
				} else if(status == 6) { // j kept unchanged, i duplicate in newGuys
					newGuys.erase(newGuys.begin() + i);
				}
				if(status > 0) {
					change = true;
					break;
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
				int nKnots = m->nKnotsIn(b);
				if( nKnots < m->multiplicity_) {
					removeFunc.insert(b);
					split( m->constDirection(), b, m->constParameter(), m->multiplicity_-nKnots, newFuncStp1 );
					break; // can only be split once by single meshrectangle insertion
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
				if( nKnots < m->multiplicity_ ) {
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

	// clear cache since mesh is now changed
	builtElementCache_ = false;

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
	std::vector<double> newKnot(p+2);
	std::copy(knot.begin(), knot.begin()+p+1, newKnot.begin()+1);
	newKnot[0] = new_knot;
	std::sort(newKnot.begin(), newKnot.begin() + p+2);
	if(constDir == 0) {
		b1 = new Basisfunction(newKnot.begin()  ,  (*b)[1].begin(),  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha1);
		b2 = new Basisfunction(newKnot.begin()+1,  (*b)[1].begin(),  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha2);
	} else if(constDir == 1) {
		b1 = new Basisfunction((*b)[0].begin(), newKnot.begin()   ,  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), newKnot.begin()+1 ,  (*b)[2].begin(), b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha2);
	} else { // insert in w
		b1 = new Basisfunction((*b)[0].begin(), (*b)[1].begin(),  newKnot.begin()   , b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), (*b)[1].begin(),  newKnot.begin()+1 , b->cp(), b->dim(), order_[0], order_[1], order_[2], b->w()*alpha2);
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


void LRSplineVolume::getGlobalKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v, std::vector<double> &knot_w) const {
	getGlobalUniqueKnotVector(knot_u, knot_v, knot_w);

	// add in duplicates where apropriate
	for(uint i=0; i<knot_u.size(); i++) {
		for(uint j=0; j<meshrect_.size(); j++) {
			if(meshrect_[j]->constDirection() == 0 && meshrect_[j]->constParameter() == knot_u[i]) {
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
			if(meshrect_[j]->constDirection() == 1 && meshrect_[j]->constParameter() == knot_v[i]) {
				for(int m=1; m<meshrect_[j]->multiplicity_; m++) {
					knot_v.insert(knot_v.begin()+i, knot_v[i]);
					i++;
				}
				break;
			}
		}
	}
	for(uint i=0; i<knot_w.size(); i++) {
		for(uint j=0; j<meshrect_.size(); j++) {
			if(meshrect_[j]->constDirection() == 2 && meshrect_[j]->constParameter() == knot_w[i]) {
				for(int m=1; m<meshrect_[j]->multiplicity_; m++) {
					knot_w.insert(knot_w.begin()+i, knot_w[i]);
					i++;
				}
				break;
			}
		}
	}
}

void LRSplineVolume::getGlobalUniqueKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v, std::vector<double> &knot_w) const {
	knot_u.clear();
	knot_v.clear();
	knot_w.clear();
	// create a huge list of all meshrectangle instances
	for(uint i=0; i<meshrect_.size(); i++) {
		if(meshrect_[i]->constDirection() == 0)
			knot_u.push_back(meshrect_[i]->constParameter());
		else if(meshrect_[i]->constDirection() == 1)
			knot_v.push_back(meshrect_[i]->constParameter());
		else
			knot_w.push_back(meshrect_[i]->constParameter());
	}

	// sort them and remove duplicates
	std::vector<double>::iterator it;
	std::sort(knot_u.begin(), knot_u.end());
	std::sort(knot_v.begin(), knot_v.end());
	std::sort(knot_w.begin(), knot_w.end());
	it = std::unique(knot_u.begin(), knot_u.end());
	knot_u.resize(it-knot_u.begin());
	it = std::unique(knot_v.begin(), knot_v.end());
	knot_v.resize(it-knot_v.begin());
	it = std::unique(knot_w.begin(), knot_w.end());
	knot_w.resize(it-knot_w.begin());
}

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
		int start[] = {-1,-1,-1};
		std::vector<std::vector<double> > row(3);
		std::vector<std::vector<double> > knot(3);
		for(int d=0; d<3; d++) {
			for(int i=0; i<order_[d]+1; i++)
				knot[d].push_back( (*b)[d][i] );
			row[d].push_back(1);
		}


		for(int d=0; d<3; d++) {

			double min = el->getParmin(d);
			double max = el->getParmax(d);
			while(knot[d][++start[d]] < min);
			while(true) {
				int p    = order_[d]-1;
				int newI = -1;
				double z;
				if(       knot[d].size() < (uint) start[d]+order_[d]   || knot[d][start[d]+  order_[d]-1] != min) {
					z    = min;
					newI = start[d];
				} else if(knot[d].size() < (uint) start[d]+2*order_[d] || knot[d][start[d]+2*order_[d]-1] != max ) {
					z    = max;
					newI = start[d] + order_[d];
				} else {
					break;
				}

				std::vector<double> newRow(row[d].size()+1, 0);
				for(uint k=0; k<row[d].size(); k++) {
					#define U(x) ( knot[d][x+k] )
					if(z < U(0) || z > U(p+1)) {
						newRow[k] = row[d][k];
						continue;
					}
					double alpha1 = (U(p) <=  z  ) ? 1 : double(   z    - U(0)) / (  U(p)  - U(0));
					double alpha2 = (z    <= U(1)) ? 1 : double( U(p+1) - z   ) / ( U(p+1) - U(1));
					newRow[k]   += row[d][k]*alpha1;
					newRow[k+1] += row[d][k]*alpha2;
					#undef U
				}
				knot[d].insert(knot[d].begin()+newI, z);
				row[d] = newRow;
			}

		}

		int ip = 0;
		for(int w=start[2]; w<start[2]+order_[2]; w++)
			for(int v=start[1]; v<start[1]+order_[1]; v++)
				for(int u=start[0]; u<start[0]+order_[0]; u++)
					for(int d=0; d<dim_; d++)
						controlPoints[ip++] += b->cp()[d]*row[0][u]*row[1][v]*row[2][w]*b->w();
	}
}

void LRSplineVolume::getBezierExtraction(int iEl, std::vector<double> &extractMatrix) const {
	Element *el = element_[iEl];
	int width  = order_[0]*order_[1]*order_[2];
	int height = el->nBasisFunctions();
	extractMatrix.clear();
	extractMatrix.resize(width*height);

	int rowI = 0;
	for(Basisfunction* b : el->support()) {
		int start[] = {-1,-1,-1};
		std::vector<std::vector<double> > row(3);
		std::vector<std::vector<double> > knot(3);
		for(int d=0; d<3; d++) {
			for(int i=0; i<order_[d]+1; i++)
				knot[d].push_back( (*b)[d][i] );
			row[d].push_back(1);
		}


		for(int d=0; d<3; d++) {

			double min = el->getParmin(d);
			double max = el->getParmax(d);
			while(knot[d][++start[d]] < min);
			while(true) {
				int p    = order_[d]-1;
				int newI = -1;
				double z;
				if(       knot[d].size() < (uint) start[d]+order_[d]   || knot[d][start[d]+  order_[d]-1] != min) {
					z    = min;
					newI = start[d];
				} else if(knot[d].size() < (uint) start[d]+2*order_[d] || knot[d][start[d]+2*order_[d]-1] != max ) {
					z    = max;
					newI = start[d] + order_[d];
				} else {
					break;
				}

				std::vector<double> newRow(row[d].size()+1, 0);
				for(uint k=0; k<row[d].size(); k++) {
					#define U(x) ( knot[d][x+k] )
					if(z < U(0) || z > U(p+1)) {
						newRow[k] = row[d][k];
						continue;
					}
					double alpha1 = (U(p) <=  z  ) ? 1 : double(   z    - U(0)) / (  U(p)  - U(0));
					double alpha2 = (z    <= U(1)) ? 1 : double( U(p+1) - z   ) / ( U(p+1) - U(1));
					newRow[k]   += row[d][k]*alpha1;
					newRow[k+1] += row[d][k]*alpha2;
					#undef U
				}
				knot[d].insert(knot[d].begin()+newI, z);
				row[d] = newRow;
			}

		}

		int colI = 0;
		for(int w=start[2]; w<start[2]+order_[2]; w++)
			for(int v=start[1]; v<start[1]+order_[1]; v++)
				for(int u=start[0]; u<start[0]+order_[0]; u++, colI++)
					extractMatrix[colI*height + rowI] += row[0][u]*row[1][v]*row[2][w]*b->w();
		rowI++;
	}
}

// compute a^n mod p where p is a prime
template <typename T>
T modPower(T a, T p, T n) {
	T result = 1;
	T tmp    = a;
	while(n>0) {
		if(n % 2)
			result = (result * tmp) % p;
		tmp = (tmp * tmp) % p;
		n /= 2;
	}
	return result;
}


// compute a^-1 mod p where p is a prime
template <typename T>
T modInverse(T a, T p) {
	return modPower(a, p, p-2);
}
bool LRSplineVolume::isLinearIndepByMappingMatrix(bool verbose) const {
#ifdef TIME_LRSPLINE
	PROFILE("Linear independent)");
#endif
	// try and figure out this thing by the projection matrix C

	std::vector<double> knots_u, knots_v, knots_w;
	getGlobalKnotVector(knots_u, knots_v, knots_w);
	int nmb_bas = basis_.size();
	int n1 = knots_u.size() - order_[0];
	int n2 = knots_v.size() - order_[1];
	int n3 = knots_w.size() - order_[2];
	int fullDim = n1*n2*n3;
	bool fullVerbose   = fullDim <  30 && nmb_bas <  50;
	bool sparseVerbose = fullDim < 250 && nmb_bas < 100;
	unsigned long long prime    = 0x7FFFFFFF;

	std::vector<std::vector<unsigned long long> > C;  // rational projection matrix

	// scaling factor to ensure that all knots are integers (assuming all multiplum of smallest knot span)
	double smallKnotU = DBL_MAX;
	double smallKnotV = DBL_MAX;
	double smallKnotW = DBL_MAX;
	for(uint i=0; i<knots_u.size()-1; i++)
		if(knots_u[i+1]-knots_u[i] < smallKnotU && knots_u[i+1] != knots_u[i])
			smallKnotU = knots_u[i+1]-knots_u[i];
	for(uint i=0; i<knots_v.size()-1; i++)
		if(knots_v[i+1]-knots_v[i] < smallKnotV && knots_v[i+1] != knots_v[i])
			smallKnotV = knots_v[i+1]-knots_v[i];
	for(uint i=0; i<knots_w.size()-1; i++)
		if(knots_w[i+1]-knots_w[i] < smallKnotW && knots_w[i+1] != knots_w[i])
			smallKnotW = knots_w[i+1]-knots_w[i];

	HashSet_const_iterator<Basisfunction*> cit_begin = basis_.begin();
	HashSet_const_iterator<Basisfunction*> cit_end   = basis_.end();
	HashSet_const_iterator<Basisfunction*> it;
	for(it=cit_begin; it != cit_end; ++it) {
		Basisfunction *b = *it;
		int startU, startV, startW;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());
		std::vector<double> locKnotW((*b)[2].begin(), (*b)[2].end());

		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == (*b)[0][0]) break;
		for(int j=0; j<order_[0]; j++) {
			if(knots_u[startU] == (*b)[0][j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == (*b)[1][0]) break;
		for(int j=0; j<order_[1]; j++) {
			if(knots_v[startV] == (*b)[1][j]) startV--;
			else break;
		}
		startV++;
		for(startW=knots_w.size(); startW-->0; )
			if(knots_w[startW] == (*b)[2][0]) break;
		for(int j=0; j<order_[2]; j++) {
			if(knots_w[startW] == (*b)[2][j]) startW--;
			else break;
		}
		startW++;

		std::vector<unsigned long long> rowU(1,1), rowV(1,1), rowW(1,1);
		int curU = startU+1;
		for(uint j=0; j<locKnotU.size()-1; j++, curU++) {
			if(locKnotU[j+1] != knots_u[curU]) {
				std::vector<unsigned long long> newRowU(rowU.size()+1, 0);
				for(uint k=0; k<rowU.size(); k++) {
					#define U(x) ((unsigned long long) (locKnotU[x+k]/smallKnotU + 0.5))
					unsigned long long z = (unsigned long long) (knots_u[curU] / smallKnotU + 0.5);
					int p = order_[0]-1;
					if(z < U(0) || z > U(p+1)) {
						newRowU[k] = rowU[k];
						continue;
					}
					unsigned long long alpha1 = (U(p) <=  z  ) ? 1 : ((  z    - U(0)) * modInverse( U(p)  - U(0), prime)  % prime);
					unsigned long long alpha2 = (z    <= U(1)) ? 1 : ((U(p+1) - z   ) * modInverse(U(p+1) - U(1), prime)  % prime);
					newRowU[k]   += rowU[k] * alpha1;
					newRowU[k+1] += rowU[k] * alpha2;
					newRowU[k]   %= prime;
					newRowU[k+1] %= prime;
					#undef U
				}
				locKnotU.insert(locKnotU.begin()+(j+1), knots_u[curU]);
				rowU = newRowU;
			}
		}
		int curV = startV+1;
		for(uint j=0; j<locKnotV.size()-1; j++, curV++) {
			if(locKnotV[j+1] != knots_v[curV]) {
				std::vector<unsigned long long> newRowV(rowV.size()+1, 0);
				for(uint k=0; k<rowV.size(); k++) {
					#define V(x) ((unsigned long long) (locKnotV[x+k]/smallKnotV + 0.5))
					unsigned long long z = (unsigned long long) (knots_v[curV] / smallKnotV + 0.5);
					int p = order_[1]-1;
					if(z < V(0) || z > V(p+1)) {
						newRowV[k] = rowV[k];
						continue;
					}
					unsigned long long alpha1 = (V(p) <=  z  ) ? 1 : ((   z    - V(0)) * modInverse(  V(p)  - V(0), prime) % prime);
					unsigned long long alpha2 = (z    <= V(1)) ? 1 : (( V(p+1) - z   ) * modInverse( V(p+1) - V(1), prime) % prime);
					newRowV[k]   += rowV[k] * alpha1;
					newRowV[k+1] += rowV[k] * alpha2;
					newRowV[k]   %= prime;
					newRowV[k+1] %= prime;
					#undef V
				}
				locKnotV.insert(locKnotV.begin()+(j+1), knots_v[curV]);
				rowV= newRowV;
			}
		}
		int curW = startW+1;
		for(uint j=0; j<locKnotW.size()-1; j++, curW++) {
			if(locKnotW[j+1] != knots_w[curW]) {
				std::vector<unsigned long long> newRowW(rowW.size()+1, 0);
				for(uint k=0; k<rowW.size(); k++) {
					#define W(x) ((unsigned long long) (locKnotW[x+k]/smallKnotW + 0.5))
					unsigned long long z = (unsigned long long) (knots_w[curW] / smallKnotW + 0.5);
					int p = order_[2]-1;
					if(z < W(0) || z > W(p+1)) {
						newRowW[k] = rowW[k];
						continue;
					}
					unsigned long long alpha1 = (W(p) <=  z  ) ? 1 : ((   z    - W(0)) * modInverse(  W(p)  - W(0), prime) % prime);
					unsigned long long alpha2 = (z    <= W(1)) ? 1 : (( W(p+1) - z   ) * modInverse( W(p+1) - W(1), prime) % prime);
					newRowW[k]   += rowW[k] * alpha1;
					newRowW[k+1] += rowW[k] * alpha2;
					newRowW[k]   %= prime;
					newRowW[k+1] %= prime;
					#undef W
				}
				locKnotW.insert(locKnotW.begin()+(j+1), knots_w[curW]);
				rowW= newRowW;
			}
		}
		std::vector<unsigned long long> totalRow(fullDim, 0);
		for(uint i1=0; i1<rowU.size(); i1++)
			for(uint i2=0; i2<rowV.size(); i2++)
				for(uint i3=0; i3<rowW.size(); i3++)
					totalRow[(startW+i3)*n1*n2 + (startV+i2)*n1 + (startU+i1)] = rowW[i3] * rowV[i2] % prime * rowU[i1] % prime;

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
		int maxI = -1;
		for(uint j=i; j<C.size(); j++) {
			if( C[j][i+zeroColumns] != 0) {
				maxI = j;
				break;
			}
		}
		if(maxI == -1) {
			i--;
			zeroColumns++;
			continue;
		}
		std::vector<unsigned long long> tmp = C[i];
		C[i] = C[maxI];
		C[maxI] = tmp;
		for(uint j=i+1; j<C.size(); j++) {
			if(C[j][i+zeroColumns] == 0) continue;
			unsigned long long scale =  C[j][i+zeroColumns] * modInverse(C[i][i+zeroColumns], prime) % prime;
			for(uint k=i+zeroColumns; k<C[j].size(); k++) {
				C[j][k] += prime - (C[i][k] * scale % prime);
				C[j][k] %= prime;
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

	int rank = (nmb_bas < n1*n2*n3-zeroColumns) ? nmb_bas : n1*n2*n3-zeroColumns;
	if(verbose) {
		std::cout << "Matrix size : " << nmb_bas << " x " << n1*n2*n3 << std::endl;
		std::cout << "Matrix rank : " << rank << std::endl;
	}

	return rank == nmb_bas;
}

void LRSplineVolume::getDiagonalElements(std::vector<int> &result) const  {
	result.clear();
	for(uint i=0; i<element_.size(); i++)
		if(element_[i]->getParmin(0) == element_[i]->getParmin(1) && element_[i]->getParmin(0) == element_[i]->getParmin(2))
			result.push_back(i);
}

void LRSplineVolume::getDiagonalBasisfunctions(std::vector<int> &result) const  {
	result.clear();
	for(Basisfunction *b : basis_) {
		bool isDiag = true;
		for(int j=0; j<order_[0]+1; j++)
			if(b->getknots(0)[j] != b->getknots(1)[j] || b->getknots(0)[j] != b->getknots(2)[j])
				isDiag = false;
		if(isDiag)
			result.push_back(b->getId());
	}
}

/************************************************************************************************************************//**
 * \brief Gets the basis functions for mixed finite element codes by generating different polynomial degree and/or continutity
 * \param raise_p1 polynomial degree to raise first parametric direction (possibly negative)
 * \param raise_p2 polynomial degree to raise second parametric direction (possibly negative)
 * \param raise_p3 polynomial degree to raise third parametric direction (possibly negative)
 * \param lower_k1 lower continuity by this amount in first parametric direction
 * \param lower_k2 lower continuity by this amount in second parametric direction
 * \param lower_k3 lower continuity by this amount in third parametric direction
 * \details This generates a brand new LRSplineVolume objects which has different polynomial degree and continuity, but
 *          meshlines in the same place.
 ***************************************************************************************************************************/
LRSplineVolume* LRSplineVolume::getDerivedBasis(int raise_p1, int raise_p2, int raise_p3, size_t lower_k1, size_t lower_k2, size_t lower_k3, int dim) const {
	// error test input
	if((raise_p1 < 0 && ((size_t) -raise_p1)>lower_k1) ||
	   (raise_p2 < 0 && ((size_t) -raise_p2)>lower_k2) ||
	   (raise_p3 < 0 && ((size_t) -raise_p3)>lower_k3) ){
	  std::cerr << "Error: getDerivedBasis undefined for raise_p < 0 and raise_p > lower_k" << std::endl;
		return NULL;
	}
	int p1 = order_[0] + raise_p1;
	int p2 = order_[1] + raise_p2;
	int p3 = order_[2] + raise_p3;
	std::vector<double> knotU(2*p1);
	std::vector<double> knotV(2*p2);
	std::vector<double> knotW(2*p3);
	for(int i=0; i<p1; i++)
		knotU[i] = start_[0];
	for(int i=0; i<p1; i++)
		knotU[i+p1] = end_[0];
	for(int i=0; i<p2; i++)
		knotV[i] = start_[1];
	for(int i=0; i<p2; i++)
		knotV[i+p2] = end_[1];
	for(int i=0; i<p3; i++)
		knotW[i] = start_[2];
	for(int i=0; i<p3; i++)
		knotW[i+p3] = end_[2];
	int N = dim * p1*p2*p3;
	std::vector<double> coef(N, 0.0);

	LRSplineVolume *result = new LRSplineVolume(p1, p2, p3, p1, p2, p3, knotU.begin(), knotV.begin(), knotW.begin(), coef.begin(), dim, false);

	for(MeshRectangle *m : meshrect_) {
		// skip end-lines
		if(m->constDirection()==0 && (fabs(m->constParameter()-start_[0])<DOUBLE_TOL || fabs(m->constParameter()-end_[0])<DOUBLE_TOL)) continue;
		if(m->constDirection()==1 && (fabs(m->constParameter()-start_[1])<DOUBLE_TOL || fabs(m->constParameter()-end_[1])<DOUBLE_TOL)) continue;
		if(m->constDirection()==2 && (fabs(m->constParameter()-start_[2])<DOUBLE_TOL || fabs(m->constParameter()-end_[2])<DOUBLE_TOL)) continue;

		int dk = 0;
		if(m->constDirection() == 0)
			dk = lower_k1 + raise_p1;
		else if(m->constDirection() == 1)
			dk = lower_k2 + raise_p2;
		else if(m->constDirection() == 2)
			dk = lower_k3 + raise_p3;

		MeshRectangle *rect = m->copy();
		rect->multiplicity_ += dk;
		result->insert_line(rect);
	}
	result->aPosterioriFixElements();
	return result;
}

int LRSplineVolume::getMinContinuity(int i) const {
	int p = order_[i];
	int minCont = p;
	for(auto rect : getAllMeshRectangles())
		if(rect->constDirection() == i)
			if(rect->multiplicity() != p) // skip C^{-1} lines (typically the edges)
				minCont = std::min(minCont, p - rect->multiplicity() - 1);
	return minCont;
}

int LRSplineVolume::getMaxContinuity(int i) const {
	int p = order_[i];
	int maxCont = -1;
	for(auto rect : getAllMeshRectangles())
		if(rect->constDirection() == i)
			if(rect->multiplicity() != p) // skip C^{-1} lines (typically the edges)
				maxCont = std::max(maxCont, p - rect->multiplicity() - 1);
	return maxCont;
}

void LRSplineVolume::setMaxContinuity(int dir, int maxCont) {
	int p = order_[dir];
	for(auto rect : getAllMeshRectangles()) {
		if(rect->constDirection() == dir) {
			int thisCont = p - rect->multiplicity() - 1;
			if(thisCont > maxCont) {
				MeshRectangle *newRect = rect->copy();
				newRect->multiplicity_ = p - maxCont - 1;
				insert_line(newRect);
			}
		}
	}
}

bool LRSplineVolume::setControlPointsVDSA(const LRSplineVolume* lr) {
	if(dim_ != lr->dimension())
		this->rebuildDimension(lr->dimension());

	std::vector<double> u(3), newCP(3);
	for(auto bit : this->getAllBasisfunctions()) {
		bit->getGrevilleParameter(u);
		std::vector<double>::iterator cp = (*bit).cp();
		lr->point(newCP, u[0], u[1], u[2]);
		for(int i=0; i<dim_; i++)
			cp[i] = newCP[i];
	}
	return true;
}

/************************************************************************************************************************//**
 * \brief functions inserting batch of rectangles (i.e. getDerivativeSpace, getPrimalSpace) may not do proper element-splits
 *        during refinement. This function fixes them a priori.
 ***************************************************************************************************************************/
void LRSplineVolume::aPosterioriFixElements() {
	for(uint i=0; i<element_.size(); i++) {
		for(uint j=0; j<meshrect_.size(); j++) {
			MeshRectangle *m = meshrect_[j];
			if(m->splits(element_[i])) {
				element_.push_back(element_[i]->split(m->constDirection(), m->constParameter()) );
				i=-1;
				break;
			}
		}
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
	basis_.clear();
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

