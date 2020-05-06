#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"

#include <set>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#ifdef HAS_BOOST
#include <boost/rational.hpp>
#endif
#include <cfloat>
#include <cmath>

typedef unsigned int uint;

//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;

namespace LR {

#define DOUBLE_TOL 1e-13


/************************************************************************************************************************//**
 * \brief Default constructor. Creates an empty LRSplineSurface object
 ***************************************************************************************************************************/
LRSplineSurface::LRSplineSurface() {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	order_.resize(2);
	start_.resize(2);
	end_.resize(2);
	rational_  = false;
	dim_       = 0;
	order_[0]  = 0;
	order_[1]  = 0;
	start_[0]  = 0;
	start_[1]  = 0;
	end_[0]    = 0;
	end_[1]    = 0;
	meshline_  = std::vector<Meshline*>(0);
	element_   = std::vector<Element*>(0);

	initMeta();
}

#ifdef HAS_GOTOOLS
/************************************************************************************************************************//**
 * \brief Copy constructor. Creates an LRSplineSurface representation from a Go::SplineSurface one
 * \param surf The SplineSurface to copy
 ***************************************************************************************************************************/
LRSplineSurface::LRSplineSurface(Go::SplineSurface *surf) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();
	initCore(surf->numCoefs_u(),      surf->numCoefs_v(),
	         surf->order_u(),         surf->order_v(),
	         surf->basis_u().begin(), surf->basis_v().begin(),
	         surf->ctrl_begin(), surf->dimension(), surf->rational() );
}
#endif


/************************************************************************************************************************//**
 * \brief Constructor. Creates a tensor product LRSplineSurface. Same signature as Go::SplineSurface constructor
 * \param n1 The number of basisfunctions in the first parametric direction
 * \param n2 The number of basisfunctions in the second parametric direction
 * \param order_u The order (polynomial degree + 1) in the first parametric direction
 * \param order_v The order (polynomial degree + 1) in the second parametric direction
 * \param knot_u The first knot vector consisting of n1+order_u elements
 * \param knot_v The second knot vector consisting of n2+order_v elements
 * \param coef Pointer to a list of n1*n2 control points describing the surface
 * \param dim The number of components of the control points, 2 for plane geometries, 3 for surfaces in 3D space
 * \param rational True if this is a NURBS surface. If so, the last component of the control points is treated as the weight
 ***************************************************************************************************************************/
/*
LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();
	initCore(n1,n2,order_u, order_v, knot_u, knot_v, coef, dim, rational);
}
*/

/************************************************************************************************************************//**
 * \brief Constructor. Creates a tensor product LRSplineSurface of the 2D unit square
 * \param n1 The number of basisfunctions in the first parametric direction
 * \param n2 The number of basisfunctions in the second parametric direction
 * \param order_u The order (polynomial degree + 1) in the first parametric direction
 * \param order_v The order (polynomial degree + 1) in the second parametric direction
 * \param knot_u The first knot vector consisting of n1+order_u elements
 * \param knot_v The second knot vector consisting of n2+order_v elements
 ***************************************************************************************************************************/
LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();

	std::vector<double> grevU = LRSpline::getGrevillePoints(order_u, knot_u, knot_u + n1 + order_u);
	std::vector<double> grevV = LRSpline::getGrevillePoints(order_v, knot_v, knot_v + n2 + order_v);
	std::vector<double> coef(n1*n2*2);
	int k=0;
	for(int j=0; j<n2; j++)
		for(int i=0; i<n1; i++) {
			coef[k++] = grevU[i];
			coef[k++] = grevV[j];
	}

	initCore(n1, n2, order_u, order_v, knot_u, knot_v, coef.begin(), 2, false);
}

/************************************************************************************************************************//**
 * \brief Constructor. Creates a uniform tensor product LRSplineSurface of the 2D unit square
 * \param n1 The number of basisfunctions in the first parametric direction
 * \param n2 The number of basisfunctions in the second parametric direction
 * \param order_u The order (polynomial degree + 1) in the first parametric direction
 * \param order_v The order (polynomial degree + 1) in the second parametric direction
 ***************************************************************************************************************************/
LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	initMeta();
	// sanity check input
	if(n1 < order_u ||
	   n2 < order_v) {
		std::cerr << "Error: n<p in LRSplineSurface constructor\n";
		// really ought to throw exception here, but havent the framework
		// for this up and running yet. Make it a zombie surface
		double knot[4] = {0,0,1,1};
		double cp[4] = {0,0,0,0};
		initCore(2,2,2,2,knot,knot,cp,1); // init dummy-state object
		return;
	}


	std::vector<double> knot_u = LRSpline::getUniformKnotVector(n1, order_u);
	std::vector<double> knot_v = LRSpline::getUniformKnotVector(n2, order_v);
	std::vector<double> grev_u = LRSpline::getGrevillePoints(order_u, knot_u.begin(), knot_u.end());
	std::vector<double> grev_v = LRSpline::getGrevillePoints(order_v, knot_v.begin(), knot_v.end());
	std::vector<double> coef(n1*n2*2);
	int k=0;
	for(int j=0; j<n2; j++)
		for(int i=0; i<n1; i++) {
			coef[k++] = grev_u[i];
			coef[k++] = grev_v[j];
	}

	// generate the uniform knot vector
	initCore(n1, n2, order_u, order_v, knot_u.begin(), knot_v.begin(), coef.begin(), 2, false);
}

void LRSplineSurface::initMeta() {
	maxTjoints_           = -1;
	doCloseGaps_          = true;
	maxAspectRatio_       = 2.0;
	doAspectRatioFix_     = false;
	refStrat_             = LR_FULLSPAN;
	refKnotlineMult_      = 1;
	symmetry_             = 1;
	element_red           = 0.5;
	element_green         = 0.5;
	element_blue          = 0.5;
	basis_red             = 1.0;
	basis_green           = 0.83;
	basis_blue            = 0.0;
	selected_basis_red    = 1.0;
	selected_basis_green  = 0.2;
	selected_basis_blue   = 0.05;
	builtElementCache_    = false;
}

/************************************************************************************************************************//**
 * \brief Destructor. Frees up all memory consumed by this LRSplineSurface
 * \details This deletes all Basisfunction, Element and Meshline used by this class. All pointers to these objects are
 *          rendered invalid
 ***************************************************************************************************************************/
LRSplineSurface::~LRSplineSurface() {
	for(Basisfunction* b : basis_)
		delete b;
	for(uint i=0; i<meshline_.size(); i++)
		delete meshline_[i];
	for(uint i=0; i<element_.size(); i++)
		delete element_[i];
}


/************************************************************************************************************************//**
 * \brief Takes a deep copy of the LRSplineSurface
 * \returns A deep copy, where all Basisfunction, Element and Meshline are also copied
 ***************************************************************************************************************************/
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
	returnvalue->order_[0]          = this->order_[0];
	returnvalue->order_[1]          = this->order_[1];
	returnvalue->start_[0]          = this->start_[0];
	returnvalue->start_[1]          = this->start_[1];
	returnvalue->end_[0]            = this->end_[0];
	returnvalue->end_[1]            = this->end_[1];
	returnvalue->maxTjoints_       = this->maxTjoints_;
	returnvalue->doCloseGaps_      = this->doCloseGaps_;
	returnvalue->doAspectRatioFix_ = this->doAspectRatioFix_;
	returnvalue->maxAspectRatio_   = this->maxAspectRatio_;

	return returnvalue;
}


#ifdef HAS_GOTOOLS
/************************************************************************************************************************//**
 * \brief Evaluate the surface at a point (u,v)
 * \param[out] pt The result, i.e. the parametric surface mapped to physical space
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param iEl The element index which this point is contained in. If used will speed up computational efficiency
 * \param u_from_right True if first coordinate should be evaluated in the limit from the right
 * \param v_from_right True if second coordinate should be evaluated in the limit from the right
 ***************************************************************************************************************************/
void LRSplineSurface::point(Go::Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const {
	std::vector<std::vector<double> > res;
	point(res, u, v, 0, u_from_right, v_from_right, iEl);
	pt = Go::Point(res[0].begin(), res[0].end());
	return ;

}

/************************************************************************************************************************//**
 * \brief Evaluate the surface at a point (u,v)
 * \param[out] pt The result, i.e. the parametric surface mapped to physical space
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param iEl The element index which this point is contained in. If used will speed up computational efficiency
 ***************************************************************************************************************************/
void LRSplineSurface::point(Go::Point &pt, double u, double v, int iEl) const {
	std::vector<std::vector<double> > res;
	point(res, u, v, 0, iEl);
	pt = Go::Point(res[0].begin(), res[0].end());
	return ;

}

/************************************************************************************************************************//**
 * \brief Evaluate the surface and its derivatives at a point (u,v)
 * \param[out] pts The result, i.e. the parametric surface as well as all parametric derivatives
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param derivs The number of derivatives requested
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details Pts is a vector where the first element is the values themselves, the next two is the first derivatives du and dv.
 *          The next three is the second derivatives in order d2u, dudv and d2v. The next four are in order d3u, d2udv, dvd2u and d3v.
 ***************************************************************************************************************************/
void LRSplineSurface::point(std::vector<Go::Point> &pts, double u, double v, int derivs, int iEl) const {

	std::vector<std::vector<double> > res;
	point(res, u, v, derivs, iEl);

	// clear and resize output array (optimization may consider this an outside task)
	pts.resize((derivs+1)*(derivs+2)/2);
	for(uint i=0; i<pts.size(); i++) {
		pts[i].resize(dim_);
		pts[i].setValue(0.0);
		for(int j=0; j<dim_; j++)
			pts[i][j] = res[i][j];
	}
	return;

}

Go::SplineCurve* LRSplineSurface::edgeCurve(parameterEdge edge,
                                            std::vector<Basisfunction*>& functions) const
{
	functions.clear();
	getEdgeFunctions(functions, edge);

	int dir = (edge == WEST || edge == EAST) ? 1 : 0;
	std::sort(functions.begin(), functions.end(),
	          [dir](const Basisfunction* a, const Basisfunction* b)
	          {
	            if (a->getGrevilleParameter()[dir] != b->getGrevilleParameter()[dir])
	              return a->getGrevilleParameter()[dir] < b->getGrevilleParameter()[dir];

	            size_t idx = 0;
	            for (auto& it : a->getknots(dir))
	              if (it != b->getknots(dir)[idx])
	                return it < b->getknots(dir)[idx++];

	            return false;
	          });

	std::vector<double> cpts;
	std::vector<double> knots;
	for (auto& it : functions) {
		std::vector<double> pt;
		it->getControlPoint(pt);
		cpts.insert(cpts.end(), pt.begin(), pt.end());
		if (knots.empty())
			knots.insert(knots.end(),it->getknots(dir).begin(),it->getknots(dir).end());
		else
			knots.push_back(it->getknots(dir).back());
	}

	return new Go::SplineCurve(Go::BsplineBasis(order(dir), knots.begin(), knots.end()),
	                           cpts.begin(), functions.front()->dim(), false);
}
#endif

/************************************************************************************************************************//**
 * \brief Evaluate the surface at a point (u,v)
 * \param[out] pt The result, i.e. the parametric surface mapped to physical space
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param iEl The element index which this point is contained in. If used will speed up computational efficiency
 * \param u_from_right True if first coordinate should be evaluated in the limit from the right
 * \param v_from_right True if second coordinate should be evaluated in the limit from the right
 ***************************************************************************************************************************/
void LRSplineSurface::point(std::vector<double> &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const {
	std::vector<std::vector<double> > res;
	point(res, u, v, 0, u_from_right, v_from_right, iEl);
	pt = res[0];
}

/************************************************************************************************************************//**
 * \brief Evaluate the surface at a point (u,v)
 * \param[out] pt The result, i.e. the parametric surface mapped to physical space
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param iEl The element index which this point is contained in. If used will speed up computational efficiency
 ***************************************************************************************************************************/
void LRSplineSurface::point(std::vector<double> &pt, double u, double v, int iEl) const {
	std::vector<std::vector<double> > res;
	point(res, u, v, 0, iEl);
	pt = res[0];
}

/************************************************************************************************************************//**
 * \brief Evaluate the surface and its derivatives at a point (u,v)
 * \param[out] pts The result, i.e. the parametric surface as well as all parametric derivatives
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param derivs The number of derivatives requested
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details Pts is a vector where the first element is the values themselves, the next two is the first derivatives du and dv.
 *          The next three is the second derivatives in order d2u, dudv and d2v. The next four are in order d3u, d2udv, dvd2u and d3v.
 ***************************************************************************************************************************/
void LRSplineSurface::point(std::vector<std::vector<double> > &pts, double u, double v, int derivs, int iEl) const {
	point(pts, u, v, derivs, u!=end_[0], v!=end_[1], iEl);
}

/************************************************************************************************************************//**
 * \brief Evaluate the surface and its derivatives at a point (u,v)
 * \param[out] pts The result, i.e. the parametric surface as well as all parametric derivatives
 * \param u The u-coordinate on which to evaluate the surface
 * \param v The v-coordinate on which to evaluate the surface
 * \param derivs The number of derivatives requested
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details Pts is a vector where the first element is the values themselves, the next two is the first derivatives du and dv.
 *          The next three is the second derivatives in order d2u, dudv and d2v. The next four are in order d3u, d2udv, dvd2u and d3v.
 ***************************************************************************************************************************/
void LRSplineSurface::point(std::vector<std::vector<double> > &pts, double u, double v, int derivs, bool u_from_right, bool v_from_right, int iEl) const {
#ifdef TIME_LRSPLINE
	PROFILE("Point()");
#endif
	std::vector<double> basis_ev;

	// clear and resize output array (optimization may consider this an outside task)
	pts.clear();
	pts.resize((derivs+1)*(derivs+2)/2);
	for(uint i=0; i<pts.size(); i++)
		pts[i].resize(dim_, 0);
	if(u < start_[0] || end_[0] < u ||
	   v < start_[1] || end_[1] < v)
		return;

	if(iEl == -1)
		iEl = getElementContaining(u,v);
	if(iEl == -1)
		return;
	for(Basisfunction* b : element_[iEl]->support() ) {
		b->evaluate(basis_ev, u,v, derivs, u_from_right, v_from_right);
		for(uint i=0; i<pts.size(); i++)
			for(int j=0; j<dim_; j++)
				pts[i][j] += basis_ev[i]*b->cp(j);
	}
}

#ifdef HAS_GOTOOLS
/************************************************************************************************************************//**
 * \brief Compute all basis functions at a parametric point (u,v)
 * \param param_u The u-coordinate on which to evaluate the surface
 * \param param_v The v-coordinate on which to evaluate the surface
 * \param[out] result All basis functions as well as parametric derivatives up to second order
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details If iEl is used, the result will contain the values of all Basisfunction active on that element (with nonzero values).
 *          If iEl is not used, then result will contain all Basisfunction in the entire LRSplineSurface object, with easier
 *          access through indexing, but will contain a lot of zeros.
 ***************************************************************************************************************************/
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
		(*it)->evaluate(values, param_u, param_v, 2, param_u!=end_[0], param_v!=end_[1]);

		result.basisValues[i]    = values[0];
		result.basisDerivs_u[i]  = values[1];
		result.basisDerivs_v[i]  = values[2];
		result.basisDerivs_uu[i] = values[3];
		result.basisDerivs_uv[i] = values[4];
		result.basisDerivs_vv[i] = values[5];
	}
}


void LRSplineSurface::computeBasis (double param_u, double param_v, Go::BasisDerivsSf3 & result, int iEl ) const {
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
		(*it)->evaluate(values, param_u, param_v, 3, param_u!=end_[0], param_v!=end_[1]);

		result.basisValues[i]    = values[0];
		result.basisDerivs_u[i]  = values[1];
		result.basisDerivs_v[i]  = values[2];
		result.basisDerivs_uu[i] = values[3];
		result.basisDerivs_uv[i] = values[4];
		result.basisDerivs_vv[i] = values[5];
		result.basisDerivs_uuu[i] = values[6];
		result.basisDerivs_uuv[i] = values[7];
		result.basisDerivs_uvv[i] = values[8];
		result.basisDerivs_vvv[i] = values[9];
	}
}

/************************************************************************************************************************//**
 * \brief Compute all basis functions at a parametric point (u,v)
 * \param param_u The u-coordinate on which to evaluate the surface
 * \param param_v The v-coordinate on which to evaluate the surface
 * \param[out] result All basis functions as well as parametric derivatives up to first order
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details If iEl is used, the result will contain the values of all Basisfunction active on that element (with nonzero values).
 *          If iEl is not used, then result will contain all Basisfunction in the entire LRSplineSurface object, with easier
 *          access through indexing, but will contain a lot of zeros.
 ***************************************************************************************************************************/
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
		(*it)->evaluate(values, param_u, param_v, 1, param_u!=end_[0], param_v!=end_[1]);

		result.basisValues[i]   = values[0];
		result.basisDerivs_u[i] = values[1];
		result.basisDerivs_v[i] = values[2];
	}
}

/************************************************************************************************************************//**
 * \brief Compute all basis functions at a parametric point (u,v)
 * \param param_u The u-coordinate on which to evaluate the surface
 * \param param_v The v-coordinate on which to evaluate the surface
 * \param[out] result All basis functions evaluated at this particular parametric point
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details If iEl is used, the result will contain the values of all Basisfunction active on that element (with nonzero values).
 *          If iEl is not used, then result will contain all Basisfunction in the entire LRSplineSurface object, with easier
 *          access through indexing, but will contain a lot of zeros.
 ***************************************************************************************************************************/
void LRSplineSurface::computeBasis(double param_u, double param_v, Go::BasisPtsSf & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	HashSet_const_iterator<Basisfunction*> it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->constSupportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->constSupportEnd();
	int nPts= (iEl<0) ? basis_.size()  : element_[iEl]->nBasisFunctions();

	result.preparePts(param_u, param_v, 0, -1, nPts);
	int i=0;
	for(it=itStart; it!=itStop; ++it, ++i)
		result.basisValues[i] = (*it)->evaluate(param_u, param_v, param_u!=end_[0], param_v!=end_[1]);
}
#endif

/************************************************************************************************************************//**
 * \brief Compute all basis functions as well as an arbitrary number of derivatives at the parametric point (u,v)
 * \param param_u The u-coordinate on which to evaluate the surface
 * \param param_v The v-coordinate on which to evaluate the surface
 * \param[out] result All basis functions and derivatives of all Basisfunction
 * \param derivs The number of derivatives requested
 * \param iEl The element index which this point is contained in. If used it will speed up computational efficiency
 * \details If iEl is used, the result will contain the values of all Basisfunction active on that element (with nonzero values).
 *          If iEl is not used, then result will contain all Basisfunction in the entire LRSplineSurface object, with easier
 *          access through indexing, but will contain a lot of zeros. The derivatives are sorted in the following 1, du, dv,
 *          d2u, dudv, d2v, d3u, d2udv, dud2v, d3v, ...
 ***************************************************************************************************************************/
void LRSplineSurface::computeBasis (double param_u,
                                    double param_v,
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
	    (*it)->evaluate(result[i], param_u, param_v, derivs, param_u!=end_[0], param_v!=end_[1]);

}

void LRSplineSurface::generateIDs() const
{
  this->LRSpline::generateIDs();
  createElementCache();
}

/************************************************************************************************************************//**
 * \brief Computes a cached lookup table for quick determination of element distribution. Allows getElementContaining()
 *        to be ran in O(log(n)) time.
 * \details This is done by creating a full tensor mesh (of elements) and storing the LR-elements from which each
 *          sub-element came from.
 ***************************************************************************************************************************/
void LRSplineSurface::createElementCache() const {
	// find the set of all unique knots in each direction (this is our global mesh)
	glob_knot_u_.clear();
	glob_knot_v_.clear();
	this->getGlobalUniqueKnotVector(glob_knot_u_, glob_knot_v_);

	// create a tensor-mesh of elements given by an nxm matrix
	elementCache_ = std::vector<std::vector<int> >(glob_knot_u_.size(), std::vector<int>(glob_knot_v_.size(), -1));

	// store parent element (i.e. the actual LR element) for each sub-cell
	for(Element *e : getAllElements()) {
		// binary-search to look up indices
		int i0 = std::lower_bound(glob_knot_u_.begin(), glob_knot_u_.end(), e->umin()) - glob_knot_u_.begin();
		int i1 = std::lower_bound(glob_knot_u_.begin(), glob_knot_u_.end(), e->umax()) - glob_knot_u_.begin();
		int j0 = std::lower_bound(glob_knot_v_.begin(), glob_knot_v_.end(), e->vmin()) - glob_knot_v_.begin();
		int j1 = std::lower_bound(glob_knot_v_.begin(), glob_knot_v_.end(), e->vmax()) - glob_knot_v_.begin();
		for(int i=i0; i<i1; i++)
			for(int j=j0; j<j1; j++)
				elementCache_[i][j] = e->getId();
	}
	builtElementCache_ = true;
}

/************************************************************************************************************************//**
 * \brief Get the element index of the element containing the parametric point (u,v)
 * \param u The u-coordinate
 * \param v The v-coordinate
 * \return The index of the element which contains (u,v)
 ***************************************************************************************************************************/
int LRSplineSurface::getElementContaining(double u, double v) const {
	// sanity check input
	if(u < startparam(0) || u > endparam(0) || v < startparam(1) || v > endparam(1))
		return -1;
	// build cache if not already present
	if(builtElementCache_ == false)
		createElementCache();

	// binary search for the right element
	size_t i = std::upper_bound(glob_knot_u_.begin(), glob_knot_u_.end(), u) - glob_knot_u_.begin() - 1;
	size_t j = std::upper_bound(glob_knot_v_.begin(), glob_knot_v_.end(), v) - glob_knot_v_.begin() - 1;

	// special case endpoints since element definitions are defined to be [umin,umax) for all except last element: [umin,umax]
	if(i == glob_knot_u_.size()-1) i--;
	if(j == glob_knot_v_.size()-1) j--;
	return elementCache_[i][j];
}

/************************************************************************************************************************//**
 * \brief Used in refinement, get the minimum span meshlines that will split at least one Basisfunction on this element
 * \param iEl The element to refine
 * \param[out] lines If the width/height ratio of the element is large, one Meshline is added to the list, if not then two
 *                   lines will be added to the list
 * \details The minimum span picks the Basisfunction with the smallest parametric support and uses this as a length measure.
 *          If several Basisfunction have equally small support, it will pick the one centered the most around the given element
 ***************************************************************************************************************************/
void LRSplineSurface::getMinspanLines(int iEl, std::vector<Meshline*>& lines) {
	Element *e = element_[iEl];
	double umin = e->umin();
	double umax = e->umax();
	double vmin = e->vmin();
	double vmax = e->vmax();
	double min_du = DBL_MAX;
	double min_dv = DBL_MAX;
	int    best_startI = order_[0]+2;
	int    best_stopI  = order_[0]+2;
	int    best_startJ = order_[1]+2;
	int    best_stopJ  = order_[1]+2;
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
	}

	if(!only_insert_span_v_line)
		lines.push_back(new Meshline(true, (e->vmin() + e->vmax())/2.0, umin, umax, 1));

	if(!only_insert_span_u_line)
		lines.push_back(new Meshline(false, (e->umin() + e->umax())/2.0, vmin, vmax, 1));

}

/************************************************************************************************************************//**
 * \brief Used in refinement, get the full span meshlines that will split all Basisfunction on this element
 * \param iEl The element to refine
 * \param[out] lines If the width/height ratio of the element is large, one Meshline is added to the list, if not then two
 *                   lines will be added to the list
 * \details The fullspan will iterate over all Basisfunction with support on this element, and use the union of their support
 *          when deciding the Meshline length.
 ***************************************************************************************************************************/
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

/************************************************************************************************************************//**
 * \brief Used in refinement, get the structured mesh meshlines that will split all Element of one Basisfunction
 * \param iBasis The Basisfunction to refine
 * \param[out] lines An unknown number of lines will be added to the list.
 * \details The structured mesh will find the largest knot span (say W) given in the local knot vector, and use W/2 as the
 *          new element sizes. It will insert new lines such that no element is larger than W/2. If any elements or local knot
 *          spans are smaller than W/2, these remain untouched.
 ***************************************************************************************************************************/
void LRSplineSurface::getStructMeshLines(Basisfunction *b, std::vector<Meshline*>& lines) {
	double umin = b->getParmin(0);
	double umax = b->getParmax(0);
	double vmin = b->getParmin(1);
	double vmax = b->getParmax(1);

	// find the largest knotspan in this function
	double max_du = 0;
	double max_dv = 0;
	for(int j=0; j<order_[0]; j++) {
		double du = (*b)[0][j+1]-(*b)[0][j];
		bool isZeroSpan =  fabs(du) < DOUBLE_TOL ;
		max_du = (isZeroSpan || max_du>du) ? max_du : du;
	}
	for(int j=0; j<order_[1]; j++) {
		double dv = (*b)[1][j+1]-(*b)[1][j];
		bool isZeroSpan =  fabs(dv) < DOUBLE_TOL ;
		max_dv = (isZeroSpan || max_dv>dv) ? max_dv : dv;
	}

	// to keep as "square" basis function as possible, only insert
	// into the largest knot spans
	for(int j=0; j<order_[0]; j++) {
		double du = (*b)[0][j+1]-(*b)[0][j];
		if( fabs(du-max_du) < DOUBLE_TOL )
			lines.push_back(new Meshline(false, ((*b)[0][j] + (*b)[0][j+1])/2.0, vmin, vmax,1));
	}
	for(int j=0; j<order_[1]; j++) {
		double dv = (*b)[1][j+1]-(*b)[1][j];
		if( fabs(dv-max_dv) < DOUBLE_TOL )
			lines.push_back(new Meshline(true, ((*b)[1][j] + (*b)[1][j+1])/2.0, umin, umax,1));
	}
}

/************************************************************************************************************************//**
 * \brief Refine one Basisfunction
 * \param index The Basisfunction to refine
 * \details This will refine a basisfunction in accordance with the structured mesh rules
 ***************************************************************************************************************************/
void LRSplineSurface::refineBasisFunction(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineBasisFunction(tmp);
}

/************************************************************************************************************************//**
 * \brief Refine several Basisfunction
 * \param indices The Basisfunction to refine
 * \details This will refine a basisfunction in accordance with the structured mesh rules. Note that the length of all lines
 *          will be precomputed before any of them are inserted. This means that you will get a different result from calling
 *          this function, rather than calling refineBasisFunction(int) several times.
 ***************************************************************************************************************************/
void LRSplineSurface::refineBasisFunction(const std::vector<int> &indices) {
	std::vector<int> sortedInd(indices);
	std::sort(sortedInd.begin(), sortedInd.end());
	std::vector<Meshline*> newLines;

	/* first retrieve all meshlines needed */
	int ib = 0;
	HashSet_iterator<Basisfunction*> it = basis_.begin();
	for(uint i=0; i<sortedInd.size(); i++) {
		while(ib < sortedInd[i]) {
			++ib;
			++it;
		}
		getStructMeshLines(*it,newLines);
	}

	/* Do the actual refinement */
	for(uint i=0; i<newLines.size(); i++) {
		Meshline *m = newLines[i];
		insert_line(!m->is_spanning_u(), m->const_par_, m->start_, m->stop_, refKnotlineMult_);
	}

	/* do a posteriori fixes to ensure a proper mesh */
	aPosterioriFixes();

	/* exit cleanly by deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++)
		delete newLines[i];
}

/************************************************************************************************************************//**
 * \brief Refine one Element
 * \param index The Element to refine
 * \details This will refine an Element in accordance with the strategy set by setRefStrat() (i.e. either fullspan or minspan)
 ***************************************************************************************************************************/
void LRSplineSurface::refineElement(int index) {
	std::vector<int> tmp = std::vector<int>(1, index);
	refineElement(tmp);
}

/************************************************************************************************************************//**
 * \brief Refine several Elements
 * \param indices The Basisfunction to refine
 * \details This will refine some Elements in accordance with the strategy set by setRefStrat() (i.e. either fullspan or minspan)
 *          Note that the length of all lines will be precomputed before any of them are inserted. This means that you will get
 *          a different result from calling this function, rather than calling refineElement(int) several times.
 ***************************************************************************************************************************/
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

	/* exit cleanly by deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++)
		delete newLines[i];
}

void LRSplineSurface::refineByDimensionIncrease(const std::vector<double> &errPerElement, double beta) {
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

	/* first retrieve all possible meshlines needed */
	std::vector<std::vector<Meshline*> > newLines(errors.size(), std::vector<Meshline*>(0));
	for(uint i=0; i<errors.size(); i++) {
		if(refStrat_ == LR_MINSPAN)
			getMinspanLines(errors[i].second, newLines[i]);
		else if(refStrat_ == LR_FULLSPAN)
			getFullspanLines(errors[i].second, newLines[i]);
		else if(refStrat_ == LR_STRUCTURED_MESH) {
			Basisfunction *b = getBasisfunction(errors[i].second);
			getStructMeshLines(b, newLines[i]);
		}
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

	/* exit cleanly by deleting all temporary new lines */
	for(uint i=0; i<newLines.size(); i++)
		for(uint j=0; j<newLines[i].size(); j++)
			delete newLines[i][j];
}

/************************************************************************************************************************//**
 * \brief Fetches list of knots that appear on a given parametric edge of the patch (used for merging multipatch models).
 * \param edge       which of the 4 parameter edges: NORTH, SOUTH, EAST or WEST
 * \param normalized set to true if the returned knots should lie in the range (0,1)
 * \returns          sorted list of mesh-lines that end on this edge
 ***************************************************************************************************************************/
std::vector<double> LRSplineSurface::getEdgeKnots(parameterEdge edge, bool normalized) const {
	std::vector<Element*> edgEl;
	this->getEdgeElements(edgEl, edge);
	std::vector<double> results;
	int pardim = edge==WEST || edge==EAST; // 0 for running u-edge, 1 for running v-edge
	for(auto el : edgEl)
		results.push_back(el->getParmin( pardim ) );
	results.push_back(this->endparam( pardim ));
	std::sort(results.begin(), results.end());
	if(normalized) {
		double u0 = this->startparam( pardim );
		double u1 = this->endparam(   pardim );
		for(uint i=0; i<results.size(); i++) {
			results[i] = (results[i] - u0) / (u1-u0);
		}
	}
	return results;
}

/************************************************************************************************************************//**
 * \brief Creates matching discretization (mesh) across patch boundaries so they can be stiched together
 * \param edge      which of the 4 parameter edges (NORTH/SOUTH/EAST/WEST) that should be refined
 * \param functions list of normalized functions (knot vectors in range (0,1)) that should appear on the boundary edge of this patch
 * \details This will refine this patche so it will have conforming mesh on the interface between this patch and another. It will
 *          not manipulate the control-points of the basis functions, and it is left to the user to make sure that these coincide
 *          if the mesh is to match in the physical space, and not just in the parametric space.
 ***************************************************************************************************************************/
void LRSplineSurface::matchParametricEdge(parameterEdge edge, const std::vector<Basisfunction*> &functions) {
	double u0 = this->startparam(0);
	double u1 = this->endparam(0);
	double v0 = this->startparam(1);
	double v1 = this->endparam(1);
	for(auto b : functions) {
		for(int d=0; d<2; d++) {
			int mult = 1;
			auto knots = b->getknots(d);
			for(uint i=0; i<knots.size(); i++) {
				if( i==knots.size()-1 || fabs(knots[i+1] - knots[i])>DOUBLE_TOL) {
					if(d==0) this->insert_line(d==0, (u1-u0)*knots[i]+u0, (v1-v0)*b->getParmin(1-d)+v0, (v1-v0)*b->getParmax(1-d)+v0, mult);
					else     this->insert_line(d==0, (v1-v0)*knots[i]+v0, (u1-u0)*b->getParmin(1-d)+u0, (u1-u0)*b->getParmax(1-d)+u0, mult);
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
bool LRSplineSurface::matchParametricEdge(parameterEdge edge, LRSplineSurface *other, parameterEdge otherEdge, bool reverse) {
	int n1 = this->nBasisFunctions();
	int n2 = other->nBasisFunctions();
	std::vector<Basisfunction*> fun1, fun2, normalized1, normalized2;
	this->getEdgeFunctions(fun1, edge);
	other->getEdgeFunctions(fun2, otherEdge);

	bool flip           = (edge==NORTH || edge==SOUTH) != (otherEdge==NORTH || otherEdge==SOUTH);
	bool reverse_along  = reverse;
	bool reverse_across = (edge==NORTH || edge==EAST)  != (otherEdge==NORTH || otherEdge==EAST);

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
	for(auto b : fun1) {
		normalized1.push_back(b->copy());
		normalized1.back()->normalize(0, this->startparam(0), this->endparam(0));
		normalized1.back()->normalize(1, this->startparam(1), this->endparam(1));
		int along  = (edge==EAST || edge==WEST);
		int across = 1-along;
		if(reverse_across)
			normalized1.back()->reverse(across);
		if(reverse_along)
			normalized1.back()->reverse(along);
		if(flip)
			normalized1.back()->flip();
	}

	// make a pre-processed rescaled set of basisfunctions that will be passed from the other patch and onto this
	for(auto b : fun2) {
		normalized2.push_back(b->copy());
		normalized2.back()->normalize(0, other->startparam(0), other->endparam(0));
		normalized2.back()->normalize(1, other->startparam(1), other->endparam(1));
		int along  = (otherEdge==EAST || otherEdge==WEST);
		int across = 1-along;
		if(reverse_across)
			normalized2.back()->reverse(across);
		if(reverse_along)
			normalized2.back()->reverse(along);
		if(flip)
			normalized2.back()->flip();
	}

	// refine patch #1
	this->matchParametricEdge(edge, normalized2);

	// refine patch #2
	other->matchParametricEdge(otherEdge, normalized1);

	// check if any change has happened:
	return (this->nBasisFunctions() != n1 || other->nBasisFunctions() != n2);

}

/************************************************************************************************************************//**
 * \brief Creates matching discretization (mesh) across patch boundaries so they can be stiched together
 * \param edge      which of the 4 parameter edges (NORTH/SOUTH/EAST/WEST) that should be refined
 * \param knots     list of normalized knots (range (0,1)) that should appear on the  boundary edge of this patch
 * \param isotropic if additional refinements should be done in an isotropic fashion (equal refinement in u- and v-direction,
 *                  i.e.function refinement)
 * \returns true if any new refinements were introduced
 * \details This will refine this patche so it will have conforming mesh on the interface between this patch and another. It will
 *          not manipulate the control-points of the basis functions, and it is left to the user to make sure that these coincide
 *          if the mesh is to match in the physical space, and not just in the parametric space.
 *          See also: getEdgeKnots()
 ***************************************************************************************************************************/
bool LRSplineSurface::matchParametricEdge(parameterEdge edge, std::vector<double> knots, bool isotropic) {
	bool didRefine = false;

	// get interface edges
	std::vector<Element*> el, el_tmp;
	this->getEdgeElements( el_tmp, edge);

	// since refinements will change or delete elements, we make a copy of the element lists
	for(auto e : el_tmp)
		el.push_back(e->copy());

	// sort elements
	int pardir = (edge==SOUTH  || edge==NORTH) ? 0 : 1;
	std::sort(el.begin(), el.end(), [pardir](const Element* a, const Element* b) {return a->getParmin(pardir) < b->getParmin(pardir);});

	// sort requested knots
	std::sort(knots.begin(), knots.end());

	// loop over elements and insert lines where needed
	int i=0;
	double u0 = this->startparam(pardir);
	double u1 = this->endparam(pardir);
	for(auto e : el) {
		double u = knots[i]*(u1-u0) + u0;
		while(fabs(u - e->getParmin(pardir) ) < DOUBLE_TOL) {
			i++;
			u = knots[i]*(u1-u0) + u0;
		}
		while(u < e->getParmax(pardir) - DOUBLE_TOL) {
			didRefine = true;
			this->insert_line(1-pardir, u, e->getParmin(1-pardir), e->getParmax(1-pardir), refKnotlineMult_);
			if(isotropic) {
				double v =  (e->getParmin(1-pardir) + e->getParmax(1-pardir)) / 2.0;
				this->insert_line(pardir, v, e->getParmin(pardir), e->getParmax(pardir), refKnotlineMult_);
			}
			i++;
			u = knots[i]*(u1-u0) + u0;
		}
	}
	return didRefine;
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
					if(fabs(left[j] - (vmin+vmax)/2) < best) {
						best = fabs(left[j] - (vmin+vmax)/2);
						bi = j;
					}
				}
				m = insert_const_v_edge(left[bi], umin, umax, refKnotlineMult_);
				if(newLines != NULL)
					newLines->push_back(m->copy());
				if(refStrat_ == LR_STRUCTURED_MESH) {
					m = insert_const_u_edge((umin+umax)/2, vmin, vmax, refKnotlineMult_);
					if(newLines != NULL)
						newLines->push_back(m->copy());
				}
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
				if(refStrat_ == LR_STRUCTURED_MESH) {
					m = insert_const_u_edge((umin+umax)/2, vmin, vmax, refKnotlineMult_);
					if(newLines != NULL)
						newLines->push_back(m->copy());
				}
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
				if(refStrat_ == LR_STRUCTURED_MESH) {
					m = insert_const_v_edge((vmin+vmax)/2, umin, umax, refKnotlineMult_);
					if(newLines != NULL)
						newLines->push_back(m->copy());
				}
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
				if(refStrat_ == LR_STRUCTURED_MESH) {
					m = insert_const_v_edge((vmin+vmax)/2, umin, umax, refKnotlineMult_);
					if(newLines != NULL)
						newLines->push_back(m->copy());
				}
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

/************************************************************************************************************************//**
 * \brief Enforces all elements to be of equal size in both directions
 * \param newLines  list of new meshlines that were inserted, or NULL if this is not needed
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
void LRSplineSurface::enforceIsotropic(std::vector<Meshline*>* newLines) {
	bool somethingFixed = true;
	while(somethingFixed) {
		somethingFixed = false;
		for(uint i=0; i<element_.size(); i++) {
			double umin = element_[i]->umin();
			double umax = element_[i]->umax();
			double vmin = element_[i]->vmin();
			double vmax = element_[i]->vmax();
			double du = umax-umin;
			double dv = vmax-vmin;
			Meshline *m;
			if(fabs(du-dv)>DOUBLE_TOL) {
				if(du>dv) {
					m = insert_line(true,  umin + du/2, vmin, vmax, refKnotlineMult_);
				} else {
					m = insert_line(false, vmin + dv/2, umin, umax, refKnotlineMult_);
				}
				if(newLines != NULL)
					newLines->push_back(m->copy());
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
	// error test input
	if(multiplicity < 1) {
		std::cerr << "LRSplineSurface::insert_line() requested line with multiplicity " << multiplicity << ". Needs non-negative values\n";
		exit(99822173);
	}
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

		if(meshline_[i]->is_spanning_u() != const_u && fabs(meshline_[i]->const_par_-const_par)<DOUBLE_TOL &&
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
			int nKnots = newline->nKnotsIn(b);
			if( nKnots < newline->multiplicity_ ) {
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
				if( nKnots < m->multiplicity_ ) {
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
	int     p                = (insert_in_u) ? order_[0] : order_[1];
	int     insert_index = 0;
	if(new_knot < knot[0] || knot[p] < new_knot)
		return ;
	while(knot[insert_index] < new_knot)
		insert_index++;
	double alpha1 = (insert_index == p)  ? 1.0 : (new_knot-knot[0])/(knot[p-1]-knot[0]);
	double alpha2 = (insert_index == 1 ) ? 1.0 : (knot[p]-new_knot)/(knot[p]-knot[1]);
	std::vector<double> newKnot(p+2);
	std::copy(knot.begin(), knot.end(), newKnot.begin() + 1);
	newKnot[0] = new_knot;
	std::sort(newKnot.begin(), newKnot.begin() + p + 2);
	if(insert_in_u) {
		b1 = new Basisfunction(newKnot.begin()  , (*b)[1].begin(), b->cp(), b->dim(), order_[0], order_[1], b->w()*alpha1);
		b2 = new Basisfunction(newKnot.begin()+1, (*b)[1].begin(), b->cp(), b->dim(), order_[0], order_[1], b->w()*alpha2);
	} else { // insert in v
		b1 = new Basisfunction((*b)[0].begin(), newKnot.begin(),     b->cp(), b->dim(), order_[0], order_[1], b->w()*alpha1);
		b2 = new Basisfunction((*b)[0].begin(), newKnot.begin() + 1, b->cp(), b->dim(), order_[0], order_[1], b->w()*alpha2);
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
			bool recursive_split = (multiplicity > 1) && ( ( insert_in_u && (*b1)[0][order_[0]]!=new_knot) ||
			                                               (!insert_in_u && (*b1)[1][order_[1]]!=new_knot)  );
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

bool LRSplineSurface::isLinearIndepByMappingMatrix(bool verbose) const {
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
	bool fullVerbose   = fullDim <  30 && nmb_bas <  50;
	bool sparseVerbose = fullDim < 250 && nmb_bas < 100;
	unsigned long long prime    = 0x7FFFFFFF;

	std::vector<std::vector<unsigned long long> > C;  // rational projection matrix

	// scaling factor to ensure that all knots are integers (assuming all multiplum of smallest knot span)
	double smallKnotU = DBL_MAX;
	double smallKnotV = DBL_MAX;
	for(uint i=0; i<knots_u.size()-1; i++)
		if(knots_u[i+1]-knots_u[i] < smallKnotU && knots_u[i+1] != knots_u[i])
			smallKnotU = knots_u[i+1]-knots_u[i];
	for(uint i=0; i<knots_v.size()-1; i++)
		if(knots_v[i+1]-knots_v[i] < smallKnotV && knots_v[i+1] != knots_v[i])
			smallKnotV = knots_v[i+1]-knots_v[i];

	HashSet_const_iterator<Basisfunction*> cit_begin = basis_.begin();
	HashSet_const_iterator<Basisfunction*> cit_end   = basis_.end();
	HashSet_const_iterator<Basisfunction*> it;
	for(it=cit_begin; it != cit_end; ++it) {
		Basisfunction *b = *it;
		int startU, startV;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());

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

		std::vector<unsigned long long> rowU(1,1), rowV(1,1);
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
		std::vector<unsigned long long> totalRow(fullDim, 0);
		for(uint i1=0; i1<rowU.size(); i1++)
			for(uint i2=0; i2<rowV.size(); i2++)
				totalRow[(startV+i2)*n1 + (startU+i1)] = rowV[i2] * rowU[i1] % prime;

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

	int rank = (nmb_bas < n1*n2-zeroColumns) ? nmb_bas : n1*n2-zeroColumns;
	if(verbose) {
		std::cout << "Matrix size : " << nmb_bas << " x " << n1*n2 << std::endl;
		std::cout << "Matrix rank : " << rank << std::endl;
	}

	return rank == nmb_bas;
}

#ifdef HAS_BOOST
void LRSplineSurface::getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const {
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

	int basCount = 0;
	for(Basisfunction *b : basis_)  {
		int startU, startV;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());

		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == (*(b))[0][0]) break;
		for(int j=0; j<order_[0]; j++) {
			if(knots_u[startU] == (*(b))[0][j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == (*(b))[1][0]) break;
		for(int j=0; j<order_[1]; j++) {
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
#endif

bool LRSplineSurface::isLinearIndepByFloatingPointMappingMatrix(bool verbose) const {
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

	for(Basisfunction *b : basis_)  {
		int startU, startV;
		std::vector<double> locKnotU((*b)[0].begin(), (*b)[0].end());
		std::vector<double> locKnotV((*b)[1].begin(), (*b)[1].end());

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

		std::vector<double > rowU(1,1), rowV(1,1);
		int curU = startU+1;
		for(uint j=0; j<locKnotU.size()-1; j++, curU++) {
			if(locKnotU[j+1] != knots_u[curU]) {
				std::vector<double > newRowU(rowU.size()+1, 0);
				for(uint k=0; k<rowU.size(); k++) {
					#define U(x) (locKnotU[x+k])
					double z = knots_u[curU];
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
		for(int j=0; j<order_[0]+1; j++)
			(*b)[0][j] = floor((*b)[0][j]/scale + 0.5);
		for(int j=0; j<order_[1]+1; j++)
			(*b)[1][j] = floor((*b)[1][j]/scale + 0.5);
	}

	// scale all LRSplineSurface values
	start_[0] = floor(start_[0]/scale + 0.5);
	start_[1] = floor(start_[1]/scale + 0.5);
	end_[0]   = floor(end_[0]  /scale + 0.5);
	end_[1]   = floor(end_[1]  /scale + 0.5);

	return scale;
}

/************************************************************************************************************************//**
 * \brief Gets the basis functions corresponding to the derivatives wrt u and v (control points must be set yourself)
 * \param diffU Derivative space for differentiation wrt u
 * \param diffV Derivative space for differentiation wrt v
 * \details This generates two brand new LRSplineSurface objects which has lower polynomial degree, and lower continuity, but
 *          meshlines in the same place. For d/du space, all constant-v meshlines have reduced continuity by 1 and the polynomial
 *          degree in u-direction is reduced by 1, while v-direction remains unchanged.
 ***************************************************************************************************************************/
std::vector<LRSplineSurface*> LRSplineSurface::getDerivativeSpace() const {
	int p1 = order_[0];
	int p2 = order_[1];
	std::vector<double> knotU(2*p1);
	std::vector<double> knotV(2*p2);
	for(int i=0; i<p1; i++)
		knotU[i] = start_[0];
	for(int i=0; i<p1; i++)
		knotU[i+p1] = end_[0];
	for(int i=0; i<p2; i++)
		knotV[i] = start_[1];
	for(int i=0; i<p2; i++)
		knotV[i+p2] = end_[1];
	int N1 = dim_ * (p1-1)*p2;
	int N2 = dim_ * p1*(p2-1);
	std::vector<double> coef1(N1);
	std::vector<double> coef2(N2);
	for(int i=0; i<N1; i++) coef1[i] = 0.0;
	for(int i=0; i<N2; i++) coef2[i] = 0.0;

	LRSplineSurface *diffU = new LRSplineSurface(p1-1, p2  , p1-1, p2  , knotU.begin()+1, knotV.begin(), coef1.begin(), dim_, false);
	LRSplineSurface *diffV = new LRSplineSurface(p1  , p2-1, p1  , p2-1, knotU.begin(), knotV.begin()+1, coef2.begin(), dim_, false);

	for(Meshline *m : meshline_) {
		if( m->span_u_line_ && (fabs(m->const_par_-start_[1])<DOUBLE_TOL || fabs(m->const_par_-end_[1])<DOUBLE_TOL)) continue;
		if(!m->span_u_line_ && (fabs(m->const_par_-start_[0])<DOUBLE_TOL || fabs(m->const_par_-end_[0])<DOUBLE_TOL)) continue;
		diffU->insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, m->multiplicity_);
		diffV->insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, m->multiplicity_);
	}
	diffU->aPosterioriFixElements();
	diffV->aPosterioriFixElements();
	std::vector<LRSplineSurface*> results(2);
	results[0] = diffU;
	results[1] = diffV;
	return results;
}

/************************************************************************************************************************//**
 * \brief Gets the basis functions for mixed finite element codes by generating different polynomial degree and/or continutity
 * \param raise_p1 polynomial degree to raise first parametric direction (possibly negative)
 * \param raise_p2 polynomial degree to raise first parametric direction (possibly negative)
 * \param lower_k1 lower continuity by this amount in first parametric direction
 * \param lower_k2 lower continuity by this amount in second parametric direction
 * \details This generates a brand new LRSplineSurface objects which has different polynomial degree and continuity, but
 *          meshlines in the same place.
 ***************************************************************************************************************************/
LRSplineSurface* LRSplineSurface::getDerivedBasis(int raise_p1, int raise_p2, size_t lower_k1, size_t lower_k2, int dim) const {
	// error test input
	if((raise_p1 < 0 && ((size_t) -raise_p1)>lower_k1) ||
	   (raise_p2 < 0 && ((size_t) -raise_p2)>lower_k2) ){
	  std::cerr << "Error: getDerivedBasis undefined for raise_p < 0 and raise_p > lower_k" << std::endl;
		return NULL;
	}
	int p1 = order_[0] + raise_p1;
	int p2 = order_[1] + raise_p2;
	std::vector<double> knotU(2*p1);
	std::vector<double> knotV(2*p2);
	for(int i=0; i<p1; i++)
		knotU[i] = start_[0];
	for(int i=0; i<p1; i++)
		knotU[i+p1] = end_[0];
	for(int i=0; i<p2; i++)
		knotV[i] = start_[1];
	for(int i=0; i<p2; i++)
		knotV[i+p2] = end_[1];
	int N = dim * p1*p2;
	std::vector<double> coef(N, 0.0);

	LRSplineSurface *result = new LRSplineSurface(p1, p2, p1, p2, knotU.begin(), knotV.begin(), coef.begin(), dim, false);

	for(Meshline *m : meshline_) {
		// skip end-lines
		if( m->span_u_line_ && (fabs(m->const_par_-start_[1])<DOUBLE_TOL || fabs(m->const_par_-end_[1])<DOUBLE_TOL)) continue;
		if(!m->span_u_line_ && (fabs(m->const_par_-start_[0])<DOUBLE_TOL || fabs(m->const_par_-end_[0])<DOUBLE_TOL)) continue;
		int dk = (m->span_u_line_) ? (lower_k2 + raise_p2) : (lower_k1 + raise_p1);
		result->insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, m->multiplicity_ + dk);
	}
	result->aPosterioriFixElements();
	return result;
}

/************************************************************************************************************************//**
 * \brief Gets the basis functions corresponding to the primal space derived from a dual space (control points must be set yourself)
 * \details This generates a brand new LRSplineSurface objects which has lower polynomial degree in both directions, and lower continuity, but
 *          is restricted to a minimum of C^0 continuity. The returned LRSpline may be used for the primal space, while *this corresponds to
 *          the dual space. The dual space can be thought of as the k-refinement of the primal space. It is important that the dual space is
 *          created first as this is of the highest degree and thus comprimise a stricter requirement on meshline length.
 ***************************************************************************************************************************/
LRSplineSurface* LRSplineSurface::getPrimalSpace() const {
	int p1 = order_[0]-1;
	int p2 = order_[1]-1;
	std::vector<double> knotU(2*p1);
	std::vector<double> knotV(2*p2);
	for(int i=0; i<p1; i++)
		knotU[i] = start_[0];
	for(int i=0; i<p1; i++)
		knotU[i+p1] = end_[0];
	for(int i=0; i<p2; i++)
		knotV[i] = start_[1];
	for(int i=0; i<p2; i++)
		knotV[i+p2] = end_[1];
	int N = dim_ * p1*p2;
	std::vector<double> coef(N);
	for(int i=0; i<N; i++) coef[i] = 0.0;

	LRSplineSurface *primal = new LRSplineSurface(p1, p2, p1, p2, knotU.begin(), knotV.begin(), coef.begin(), dim_, false);

	for(Meshline *m : meshline_) {
		if( m->span_u_line_ && (fabs(m->const_par_-start_[1])<DOUBLE_TOL || fabs(m->const_par_-end_[1])<DOUBLE_TOL)) continue;
		if(!m->span_u_line_ && (fabs(m->const_par_-start_[0])<DOUBLE_TOL || fabs(m->const_par_-end_[0])<DOUBLE_TOL)) continue;
		int mult = m->multiplicity_;
		if( m->span_u_line_ && mult >= p2) mult = p2-1;
		if(!m->span_u_line_ && mult >= p1) mult = p1-1;
		primal->insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, mult);
	}
	primal->aPosterioriFixElements();
	return primal;
}

/************************************************************************************************************************//**
 * \brief Gets the basis functions corresponding to an order elevation (control points must be set yourself)
 * \param raiseOrderU The number of degrees to raise the first parametric direction
 * \param raiseOrderV The number of degrees to raise the second parametric direction
 * \returns pointer to the newly created LRSplineSurface object
 * \details This generates a brand new LRSplineSurface object which has higher polynomial degree, and meshlines in the same
 *          place, but with higher multiplicity which gives them the same continuity as the initial basis.
 ***************************************************************************************************************************/
LRSplineSurface* LRSplineSurface::getRaiseOrderSpace(int raiseOrderU, int raiseOrderV) const {
	int p1 = order_[0]+raiseOrderU;
	int p2 = order_[1]+raiseOrderV;
	std::vector<double> knotU(2*p1);
	std::vector<double> knotV(2*p2);
	for(int i=0; i<p1; i++)
		knotU[i] = start_[0];
	for(int i=0; i<p1; i++)
		knotU[i+p1] = end_[0];
	for(int i=0; i<p2; i++)
		knotV[i] = start_[1];
	for(int i=0; i<p2; i++)
		knotV[i+p2] = end_[1];
	int N = dim_ * p1*p2;
	std::vector<double> coef(N);
	for(int i=0; i<N; i++)
		coef[i] = 0.0;
	LRSplineSurface *result = new LRSplineSurface(p1,p2,p1,p2,knotU.begin(), knotV.begin(), coef.begin(), dim_, false);

	for(Meshline *m : meshline_) {
		int newMult = m->multiplicity_ + ((m->span_u_line_) ? raiseOrderV : raiseOrderU);
		result->insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, newMult);
	}
	result->aPosterioriFixElements();

	return result;
}

int LRSplineSurface::getMinContinuity(int i) const {
	int p = order_[i];
	int minCont = p;
	for(auto line : getAllMeshlines())
		if(1-line->is_spanning_u() == i)
			if(line->multiplicity() != p) // skip C^{-1} lines (typically the edges)
				minCont = std::min(minCont, p - line->multiplicity() - 1);
	return minCont;
}

int LRSplineSurface::getMaxContinuity(int i) const {
	int p = order_[i];
	int maxCont = -1;
	for(auto line : getAllMeshlines())
		if(1-line->is_spanning_u() == i)
			if(line->multiplicity() != p) // skip C^{-1} lines (typically the edges)
				maxCont = std::max(maxCont, p - line->multiplicity() - 1);
	return maxCont;
}

/************************************************************************************************************************//**
 * \brief provides a maximum regularity constraint on the entire mesh
 * \param contU The maximum continuity on all meshlines with constant u
 * \param contV The maximum continuity on all meshlines with constant v
 * \returns True if operation was successful. False if error on input parameters
 * \details Ensures that all meshlines have at least multiplictity p-c, where p is the polynomial degree and c is the requested
 *          continuity. Note that any existing meshlines of lower continuity is unchanged by this call.
 ***************************************************************************************************************************/
bool LRSplineSurface::setGlobalContinuity(int contU, int contV) {
	if(contU < -1 || contV < -1)
		return false;
	std::vector<Meshline*> existingLines;
	for(Meshline *m : meshline_)
		existingLines.push_back(m->copy());

	for(Meshline *m : existingLines) {
		int newMult = (m->span_u_line_) ? order_[1]-contV-1 : order_[0]-contU-1;
		if(newMult < 1) continue;
		insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, newMult);
	}

	// clean up
	for(Meshline *m : existingLines)
		delete m;
	return true;
}

/************************************************************************************************************************//**
 * \brief Decreases the continuity of all meshlines by fixed amount down to minimum of C^{-1}
 * \param du The amount of continuity reduction on all lines with constant u parameter
 * \param dv The amount of continuity reduction on all lines with constant v parameter
 * \returns True if operation was successful. False if error on input parameters
 ***************************************************************************************************************************/
bool LRSplineSurface::decreaseContinuity(int du, int dv) {
	if(du < 0 || dv < 0) {
		return false;
	}
	std::vector<Meshline*> existingLines;
	for(Meshline *m : meshline_)
		existingLines.push_back(m->copy());

	for(Meshline *m : existingLines) {
		int newMult = m->multiplicity_ + ((m->span_u_line_) ? dv : du);
		int maxMult = ((m->span_u_line_) ? order_[1] : order_[0]);
		if(newMult > maxMult)
			newMult = maxMult;
		insert_line(!m->span_u_line_, m->const_par_, m->start_, m->stop_, newMult);
	}

	// clean up
	for(Meshline *m : existingLines)
		delete m;
	return true;
}

bool LRSplineSurface::setControlPointsVDSA(const LRSplineSurface* lr) {
	if(dim_ != lr->dimension())
		this->rebuildDimension(lr->dimension());

	std::vector<double> u(2), newCP(2);
	for(auto bit : this->getAllBasisfunctions()) {
		bit->getGrevilleParameter(u);
		std::vector<double>::iterator cp = (*bit).cp();
		lr->point(newCP, u[0], u[1]);
		for(int i=0; i<dim_; i++)
			cp[i] = newCP[i];
	}
	return true;
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
		for(int j=0; j<order_[0]+1; j++)
			if((*b)[0][j] != (*b)[1][j])
				isDiag = false;
		if(isDiag)
			result.push_back(i);
		i++;
	}
}

/************************************************************************************************************************//**
 * \brief functions inserting batch of lines (i.e. getDerivativeSpace, getPrimalSpace) may not do proper element-splits
 *        during refinement. This function fixes them a priori.
 ***************************************************************************************************************************/
void LRSplineSurface::aPosterioriFixElements() {
	for(uint i=0; i<element_.size(); i++) {
		for(uint j=0; j<meshline_.size(); j++) {
			if(meshline_[j]->splits(element_[i])) {
				element_.push_back(element_[i]->split(meshline_[j]->is_spanning_u(), meshline_[j]->const_par_));
				i=-1;
				break;
			}
		}
	}
}

void LRSplineSurface::getBezierElement(int iEl, std::vector<double> &controlPoints) const {
	controlPoints.clear();
	controlPoints.resize(order_[0]*order_[1]*dim_, 0);
	Element *el = element_[iEl];
	for(Basisfunction* b : el->support()) {
		std::vector<double> knotU((*b)[0].begin(), (*b)[0].end() );
		std::vector<double> knotV((*b)[1].begin(), (*b)[1].end() );
		int startU=-1;
		int startV=-1;

		std::vector<double> rowU(1,1), rowV(1,1);

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

		int ip = 0;
		for(int v=startV; v<startV+order_[1]; v++)
			for(int u=startU; u<startU+order_[0]; u++)
				for(int d=0; d<dim_; d++)
					controlPoints[ip++] += b->cp()[d]*rowU[u]*rowV[v]*b->w();

	}
}

void LRSplineSurface::getBezierExtraction(int iEl, std::vector<double> &extractMatrix) const {
	Element *el = element_[iEl];
	int width  = order_[0]*order_[1];
	int height = el->nBasisFunctions();
	extractMatrix.clear();
	extractMatrix.resize(width*height);

	int rowI = 0;
	for(Basisfunction* b : el->support()) {
		int start[] = {-1,-1};
		std::vector<std::vector<double> > row(2);
		std::vector<std::vector<double> > knot(2);
		for(int d=0; d<2; d++) {
			for(int i=0; i<order_[d]+1; i++)
				knot[d].push_back( (*b)[d][i] );
			row[d].push_back(1);
		}


		for(int d=0; d<2; d++) {

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
		for(int v=start[1]; v<start[1]+order_[1]; v++)
			for(int u=start[0]; u<start[0]+order_[0]; u++, colI++)
				extractMatrix[colI*height + rowI] += row[0][u]*row[1][v]*b->w();
		rowI++;
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
	start_[0] =  DBL_MAX;
	end_[0]   = -DBL_MAX;
	start_[1] =  DBL_MAX;
	end_[1]   = -DBL_MAX;

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
	is >> order_[0];    ws(is);
	is >> order_[1];    ws(is);
	is >> nBasis;      ws(is);
	is >> nMeshlines;  ws(is);
	is >> nElements;   ws(is);
	is >> dim_;        ws(is);
	is >> rational_;   ws(is);

	meshline_.resize(nMeshlines);
	element_.resize(nElements);
	basis_.clear();
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
		Basisfunction *b = new Basisfunction(dim_, order_[0], order_[1]);
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
		start_[0] = (element_[i]->umin() < start_[0]) ? element_[i]->umin() : start_[0];
		end_[0]   = (element_[i]->umax() > end_[0]  ) ? element_[i]->umax() : end_[0]  ;
		start_[1] = (element_[i]->vmin() < start_[1]) ? element_[i]->vmin() : start_[1];
		end_[1]   = (element_[i]->vmax() > end_[1]  ) ? element_[i]->vmax() : end_[1]  ;
	}
}

void LRSplineSurface::write(std::ostream &os) const {
	generateIDs();
	os << std::setprecision(16);
	os << "# LRSPLINE SURFACE\n";
	os << "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n\t";
	os << order_[0] << "\t";
	os << order_[1] << "\t";
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
	int dx = end_[0] - start_[0];
	int dy = end_[1] - start_[1];
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;
	// set the duplicate-knot-line (dkl) display width
	double dkl_range = (min_span_u>min_span_v) ? min_span_v*scale/6.0 : min_span_u*scale/6.0;
	int xmin = (start_[0] - dx/30.0)*scale;
	int ymin = (start_[1] - dy/30.0)*scale;
	int xmax = (end_[0]   + dx/30.0)*scale + dkl_range;
	int ymax = (end_[1]   + dy/30.0)*scale + dkl_range;

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
				if(meshline_[i]->stop_ == end_[0])
					out << meshline_[i]->stop_*scale+dkl_range << " " << meshline_[i]->const_par_*scale + dm*m << " lineto\n";
				else
					out << meshline_[i]->stop_*scale << " " << meshline_[i]->const_par_*scale + dm*m << " lineto\n";
			} else {
				out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->start_*scale << " moveto\n";
				if(meshline_[i]->stop_ == end_[1])
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

				std::vector<double> pt;
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


		std::vector<double> pt;
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
		int textX   = ((*b)[0][1] == (*b)[0][order_[0]]) ? -2 : 1;
		int textY   = ((*b)[1][1] == (*b)[1][order_[1]]) ? -2 : 1;
		// move text on edge functions
		if((*b)[0][1] == end_[0])
			textX = 1;
		else if((*b)[0][order_[0]-1] == start_[0])
			textX = -2;
		if((*b)[1][1] == end_[1])
			textY = 1;
		else if((*b)[1][order_[1]-1] == start_[1])
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

	int dx = end_[0] - start_[0];
	int dy = end_[1] - start_[1];
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;

	double max_du = 0;
	double max_dv = 0;
	for(Basisfunction *b : basis_) {
		double du    = (*b)[0][order_[0]] - (*b)[0][0];
		double dv    = (*b)[1][order_[1]] - (*b)[1][0];
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
		double du    = (*b)[0][order_[0]] - (*b)[0][0];
		double dv    = (*b)[1][order_[1]] - (*b)[1][0];

		// move C^{-1} text on internal functions
		double textOffset = 15.0;
		int textX   = ((*b)[0][1] == (*b)[0][order_[0]]) ? -2 : 1;
		int textY   = ((*b)[1][1] == (*b)[1][order_[1]]) ? -2 : 1;
		// move text on edge functions
		if((*b)[0][1] == end_[0])
			textX = 1;
		else if((*b)[0][order_[0]-1] == start_[0])
			textX = -2;
		if((*b)[1][1] == end_[1])
			textY = 1;
		else if((*b)[1][order_[1]-1] == start_[1])
			textY = -2;

		for(int j=1; j<order_[0]; j++)
			avg_u += (*b)[0][j];
		for(int j=1; j<order_[1]; j++)
			avg_v += (*b)[1][j];
		avg_u /= (order_[0]-1);
		avg_v /= (order_[1]-1);

		bool doColor = false;
		if(colorBasis != NULL)
			for(uint j=0; j<colorBasis->size(); j++)
				if(i == (int) colorBasis->at(j))
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

#undef DOUBLE_TOL

} // end namespace LR

