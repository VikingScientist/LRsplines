#ifndef LRSPLINESURFACE_H
#define LRSPLINESURFACE_H

#include <vector>
#ifdef HAS_GOTOOLS
	#include <GoTools/utils/Point.h>
	#include <GoTools/geometry/SplineSurface.h>
	#include <GoTools/geometry/SplineCurve.h>
#endif
#ifdef HAS_BOOST
#include <boost/rational.hpp>
#endif
#ifdef TIME_LRSPLINE
#include "Profiler.h"
#endif
#include "LRSpline.h"
#include "Basisfunction.h"
#include "Meshline.h"
#include "Element.h"
#include "HashSet.h"

#ifdef HAS_GOTOOLS
#ifndef GOTOOLS_HAS_BASISDERIVS_SF3
namespace Go {
/// Structure for storage of results of grid evaluation of the basis function of a spline surface.
/// Position, first, second and third derivatives
struct BasisDerivsSf3
{
  /// Parameter double in which the basis functions are evaluated
  double param[2];
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
  int left_idx[2];
  /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
  std::vector< double > basisValues;

  /// the derivative of all basis functions in u direction, same size as previous
  std::vector< double > basisDerivs_u;
  /// the derivative of all basis functions in v direction, same size as previous
  std::vector< double > basisDerivs_v;

  /// the second derivative of all basis functions twice in u direction, same size as previous
  std::vector< double > basisDerivs_uu;
  /// the second derivative of all basis functions in u and v direction, same size as previous
  std::vector< double > basisDerivs_uv;
  /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;

    std::vector< double > basisDerivs_uuu;
    std::vector< double > basisDerivs_uuv;
    std::vector< double > basisDerivs_uvv;
    std::vector< double > basisDerivs_vvv;

    void prepareDerivs(double u, double v, int idx_u, int idx_v, int size)
    {
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_uu.resize(size);
      basisDerivs_uv.resize(size);
      basisDerivs_vv.resize(size);
      basisDerivs_uuu.resize(size);
      basisDerivs_uuv.resize(size);
      basisDerivs_uvv.resize(size);
      basisDerivs_vvv.resize(size);
    }
};
}
#endif
#endif


namespace LR {


class LRSplineSurface : public LRSpline {

public:
	// constructors and destructors
	LRSplineSurface();
#ifdef HAS_GOTOOLS
	explicit LRSplineSurface(Go::SplineSurface *surf);
#endif
	LRSplineSurface(int n1, int n2, int order_u, int order_v);
	LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v);
	// LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational=false);
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3>
	LRSplineSurface(int n1, int n2, int order_u, int order_v, RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 coef, int dim, bool rational=false) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
		initMeta();
		initCore(n1, n2, order_u, order_v, knot_u, knot_v, coef, dim, rational);
	}
	virtual ~LRSplineSurface();
	LRSplineSurface* copy() const;

	virtual void generateIDs() const;

	// surface evaluation
#ifdef HAS_GOTOOLS
	virtual void point(Go::Point &pt, double u, double v, int iEl=-1) const;
	virtual void point(Go::Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const;
	virtual void point(std::vector<Go::Point> &pts, double upar, double vpar, int derivs, int iEl=-1) const;
	void computeBasis (double param_u, double param_v, Go::BasisPtsSf     & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf  & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf2 & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf3 & result, int iEl=-1 ) const;

	Go::SplineCurve* edgeCurve(parameterEdge edge,
	                           std::vector<Basisfunction*>& functions) const;
#endif
	virtual void point(std::vector<double> &pt, double u, double v, int iEl=-1) const;
	virtual void point(std::vector<double> &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const;
	virtual void point(std::vector<std::vector<double> > &pts, double upar, double vpar, int derivs, int iEl=-1) const;
	virtual void point(std::vector<std::vector<double> > &pts, double upar, double vpar, int derivs, bool u_from_right, bool v_from_right, int iEl=-1) const;
	void computeBasis (double param_u,
	                   double param_v,
	                   std::vector<std::vector<double> >& result,
	                   int derivs=0,
	                   int iEl=-1 ) const;
	int getElementContaining(double u, double v) const;
	// TODO: get rid of the iEl argument in evaluation signatures - it's too easy to mess it up (especially with derivatives at multiple-knot boundaries).
	//       Try and sort the Elements after all refinements and binary search for the containing point in logarithmic time
	int getElementContaining(const std::vector<double>& parvalues) const { return this->getElementContaining(parvalues[0], parvalues[1]); };

	// refinement functions
	void refineBasisFunction(int index);
	void refineBasisFunction(const std::vector<int> &indices);
	void refineElement(int index);
	void refineElement(const std::vector<int> &indices);
	void refineByDimensionIncrease(const std::vector<double> &error, double beta);
	bool matchParametricEdge(parameterEdge edge, std::vector<double> knots, bool isotropic=false);
	void matchParametricEdge(parameterEdge edge, const std::vector<Basisfunction*> &functions);
	bool matchParametricEdge(parameterEdge edge, LRSplineSurface *other, parameterEdge otherEdge, bool reverse);

	// (private) refinement functions
	Meshline* insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity=1);
	Meshline* insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity=1);
	void getFullspanLines(  int iEl,          std::vector<Meshline*>& lines);
	void getMinspanLines(   int iEl,          std::vector<Meshline*>& lines);
	void getStructMeshLines(Basisfunction *b, std::vector<Meshline*>& lines);
	void aPosterioriFixes() ;
	void closeGaps(            std::vector<Meshline*>* newLines=NULL);
	void enforceMaxTjoints(    std::vector<Meshline*>* newLines=NULL);
	void enforceMaxAspectRatio(std::vector<Meshline*>* newLines=NULL);
	void enforceIsotropic(     std::vector<Meshline*>* newLines=NULL);

	// linear independence methods
	bool isLinearIndepByOverloading(bool verbose) ;
	bool isLinearIndepByFloatingPointMappingMatrix(bool verbose) const ;
#ifdef HAS_BOOST
	void getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const ;
#endif
	bool isLinearIndepByMappingMatrix(bool verbose) const;

	void updateSupport(Basisfunction *f) ;
	void updateSupport(Basisfunction *f,
	                   std::vector<Element*>::iterator start,
	                   std::vector<Element*>::iterator end ) ;

	// common get methods
	void getGlobalKnotVector      (std::vector<double> &knot_u, std::vector<double> &knot_v) const;
	void getGlobalUniqueKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const;
	void getBezierElement         (int iEl, std::vector<double> &controlPoints)              const;
	void getBezierExtraction      (int iEl, std::vector<double> &extractMatrix)              const;
	int nMeshlines()              const                { return meshline_.size(); };

	// more get-methods
	std::vector<Meshline*>::iterator       meshlineBegin()         { return meshline_.begin(); };
	std::vector<Meshline*>::iterator       meshlineEnd()           { return meshline_.end(); };
	Meshline* getMeshline(int i)                                   { return meshline_[i]; };
	const Meshline* getMeshline(int i) const                       { return meshline_[i]; };
	std::vector<Meshline*> getAllMeshlines()                       { return meshline_;    };
	const std::vector<Meshline*> getAllMeshlines() const           { return meshline_;    };
	int getMaxContinuity(int i) const ;
	int getMinContinuity(int i) const ;

	// assorted specialized functions
	std::vector<double> getEdgeKnots(parameterEdge edge, bool normalized=false) const;
	double makeIntegerKnots();
	void getSupportElements(       std::vector<int> &result, const std::vector<int> &basisfunctions) const ;
	void getDiagonalElements(      std::vector<int> &result) const ;
	void getDiagonalBasisfunctions(std::vector<int> &result) const ;
	void printElements(std::ostream &out) const;

	// fetch function spaces of different order/continuity
	LRSplineSurface*              getRaiseOrderSpace(int raiseOrderU, int raiseOrderV) const;
	std::vector<LRSplineSurface*> getDerivativeSpace() const ;
  LRSplineSurface*              getDerivedBasis(int raise_p1, int raise_p2, size_t lower_k1, size_t lower_k2, int dim=1) const;
	LRSplineSurface*              getPrimalSpace() const ;
	bool setGlobalContinuity(int contU, int contV);
	bool decreaseContinuity( int du,    int dv);

	//! \brief sets the controlpoints of *this to be a variational diminishing spline approximation of the input argument
	bool setControlPointsVDSA(const LRSplineSurface *lr);

	// input output methods
	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

	// print LR splines as eps-files
	void setElementColor(double r, double g, double b) ;
	void setBasisColor(double r, double g, double b) ;
	void setSelectedBasisColor(double r, double g, double b) ;
	void writePostscriptMesh(std::ostream &out, bool close=true, std::vector<int> *colorElements=NULL) const;
	void writePostscriptElements(std::ostream &out, int nu=2, int nv=2, bool close=true, std::vector<int> *colorElements=NULL) const;
	void writePostscriptFunctionSpace(std::ostream &out, std::vector<int> *colorBasis=NULL, bool drawAll=true, bool close=true) const;
	void writePostscriptMeshWithControlPoints(std::ostream &out, int nu=2, int nv=2) const ;

private:
	// make non-copyable
	LRSplineSurface(const LRSplineSurface&) = delete;
	LRSplineSurface& operator=(const LRSplineSurface&) = delete;

	// caching stuff
	mutable std::vector<std::vector<int> > elementCache_;
	mutable std::vector<double>            glob_knot_u_;
	mutable std::vector<double>            glob_knot_v_;
	mutable bool                           builtElementCache_;

	void createElementCache() const;

	// initializeation methods (called from constructors)
	void initMeta();
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3>
	void initCore(int n1, int n2, int order_u, int order_v, RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 coef, int dim, bool rational=false) {
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
		int p1     = order_u;
		int p2     = order_v;
		order_.resize(2);
		start_.resize(2);
		end_.resize(2);
		rational_  = rational;
		dim_       = dim;
		order_[0]  = order_u;
		order_[1]  = order_v;
		start_[0]  = knot_u[0];
		start_[1]  = knot_v[0];
		end_[0]    = knot_u[n1+p1-1];
		end_[1]    = knot_v[n2+p2-1];

		int unique_u=0;
		int unique_v=0;
		std::vector<int> elm_u; // element index of the knot vector (duplicate knots is multiple index in knot_u, single index in elmRows)
		std::vector<int> elm_v;
		for(int i=0; i<n1+p1; i++) {// const u, spanning v
			int mult = 1;
			elm_u.push_back(unique_u);
			while(i+1<n1+p1 && knot_u[i]==knot_u[i+1]) {
				i++;
				mult++;
				elm_u.push_back(unique_u);
			}
			unique_u++;
			meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2+p2-1], mult) );
		}
		for(int i=0; i<n2+p2; i++) {// const v, spanning u
			int mult = 1;
			elm_v.push_back(unique_v);
			while(i+1<n2+p2 && knot_v[i]==knot_v[i+1]) {
				i++;
				mult++;
				elm_v.push_back(unique_v);
			}
			unique_v++;
			meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n1+p1-1], mult) );
		}
		std::vector<std::vector<Element*> > elmRows(unique_v-1);
		for(int j=0; j<unique_v-1; j++) {
			for(int i=0; i<unique_u-1; i++) {
				double umin = meshline_[i]->const_par_;
				double vmin = meshline_[unique_u + j]->const_par_;
				double umax = meshline_[i+1]->const_par_;
				double vmax = meshline_[unique_u + j+1]->const_par_;
				Element* elm = new Element(umin, vmin, umax, vmax);
				element_.push_back(elm);
				elmRows[j].push_back(elm);
			}
		}

		for(int j=0; j<n2; j++)
			for(int i=0; i<n1; i++) {
				Basisfunction *b = new Basisfunction(knot_u+i, knot_v+j, coef+(j*n1+i)*(dim+rational), dim, order_u, order_v);
				basis_.insert(b);
				for(int k=elm_v[j]; k<elm_v[j+p2]; k++)
					updateSupport(b, elmRows[k].begin() + elm_u[i], elmRows[k].begin() + elm_u[i+p1]);
		}
	}

	void aPosterioriFixElements();
	void split(bool insert_in_u, Basisfunction* b, double new_knot, int multiplicity, HashSet<Basisfunction*> &newFunctions);
	Meshline* insert_line(bool const_u, double const_par, double start, double stop, int multiplicity);

	std::vector<Meshline*> meshline_;

	// plotting parameters
	double element_red;
	double element_green;
	double element_blue;
	double basis_red;
	double basis_green;
	double basis_blue;
	double selected_basis_red;
	double selected_basis_green;
	double selected_basis_blue;

};

} // end namespace LR

#endif
