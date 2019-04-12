#ifndef LRSPLINEVOLUME_H
#define LRSPLINEVOLUME_H

#ifdef HAS_GOTOOLS
	#include <GoTools/utils/Point.h>
	#include <GoTools/trivariate/SplineVolume.h>
#endif
#ifdef HAS_BOOST
	#include <boost/rational.hpp>
#endif
#ifdef TIME_LRSPLINE
#include "Profiler.h"
#endif
#include "LRSpline.h"
#include "HashSet.h"
#include "Basisfunction.h"
#include "MeshRectangle.h"
#include "Element.h"
#include <set>

namespace Go {
	class SplineVolume;
}

namespace LR {

class MeshRectangle;

class LRSplineVolume : public LRSpline {

public:
	LRSplineVolume();
#ifdef HAS_GOTOOLS
	LRSplineVolume(Go::SplineVolume *surf);
#endif
	LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w);
	LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, double *knot_u, double *knot_v, double *knot_w);
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3,
	          typename RandomIterator4>
	LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 knot_w, RandomIterator4 coef, int dim, bool rational = false) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
		initMeta();
		initCore(n1, n2, n3, order_u, order_v, order_w, knot_u, knot_v, knot_w, coef, dim, rational);
	}
	// LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, double *knot_u, double *knot_v, double *knot_w, double *coef, int dim, bool rational=false);
	virtual ~LRSplineVolume();

	LRSplineVolume* copy() const;
	// surface evaluation
#ifdef HAS_GOTOOLS
	virtual void point(Go::Point &pt, double u, double v, double w, int iEl=-1) const;
	virtual void point(Go::Point &pt, double u, double v, double w, int iEl, bool u_from_right, bool v_from_right, bool w_from_right) const;
	virtual void point(std::vector<Go::Point> &pts, double upar, double vpar, double wpar, int derivs, int iEl=-1) const;
	void computeBasis (double param_u, double param_v, double param_w, Go::BasisPts     & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, double param_w, Go::BasisDerivs  & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, double param_w, Go::BasisDerivs2 & result, int iEl=-1 ) const;
#endif
	virtual void point(std::vector<double> &pt, double u, double v, double w, int iEl=-1) const;
	virtual void point(std::vector<double> &pt, double u, double v, double w, int iEl, bool u_from_right, bool v_from_right, bool w_from_right) const;
	virtual void point(std::vector<std::vector<double> > &pts, double u, double v, double w, int derivs, int iEl=-1) const;
	virtual void point(std::vector<std::vector<double> > &pts, double u, double v, double w, int derivs, bool u_from_right, bool v_from_right, bool w_from_right, int iEl=-1) const;
	void computeBasis (double param_u,
	                   double param_v,
	                   double param_w,
	                   std::vector<std::vector<double> >& result,
	                   int derivs=0,
	                   int iEl=-1 ) const;
	std::set<int> getElementNeighbours(int iEl, parameterEdge edge) const;
	int getElementContaining(double u, double v, double w) const;
	// TODO: get rid of the iEl argument in evaluation signatures - it's too easy to mess it up (especially with derivatives at multiple-knot boundaries).
	//       Try and sort the Elements after all refinements and binary search for the containing point in logarithmic time
	int getElementContaining(const std::vector<double>& parvalues) const { return this->getElementContaining(parvalues[0], parvalues[1], parvalues[2]); };

	// refinement functions
	void refineElement(int index);
	void refineElement(const std::vector<int> &indices);
	void refineBasisFunction(int index);
	void refineBasisFunction(const std::vector<int> &indices);
	void refineByDimensionIncrease(const std::vector<double> &error, double beta);
	void matchParametricEdge(parameterEdge edge, const std::vector<Basisfunction*> &functions);
	bool matchParametricEdge(parameterEdge edge, LRSplineVolume *other, parameterEdge otherEdge, bool reverse_u, bool reverse_v, bool flip_uv);

	// (private) refinement functions
	void getFullspanRects(  int iEl,    std::vector<MeshRectangle*>& rects);
	void getMinspanRects(   int iEl,    std::vector<MeshRectangle*>& rects);
	void getStructMeshRects(Basisfunction *b, std::vector<MeshRectangle*>& rects);
	/*
	MeshRectangle* insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity=1);
	MeshRectangle* insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity=1);
	void aPosterioriFixes() ;
	void closeGaps(            std::vector<MeshRectangle*>* newLines=NULL);
	void enforceMaxTjoints(    std::vector<MeshRectangle*>* newLines=NULL);
	void enforceMaxAspectRatio(std::vector<MeshRectangle*>* newLines=NULL);
	*/
	bool enforceIsotropic();

	// linear independence methods
	bool isLinearIndepByOverloading(bool verbose) ;
	bool isLinearIndepByMappingMatrix(bool verbose) const ;
	bool isLinearIndepByFloatingPointMappingMatrix(bool verbose) const ;
#ifdef HAS_BOOST
	void getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const ;
#endif

	void updateSupport(Basisfunction *f) ;
	void updateSupport(Basisfunction *f,
	                   std::vector<Element*>::iterator start,
	                   std::vector<Element*>::iterator end ) ;

	// common get methods
	void getGlobalKnotVector      (std::vector<double> &knot_u,
	                               std::vector<double> &knot_v,
	                               std::vector<double> &knot_w) const;
	void getGlobalUniqueKnotVector(std::vector<double> &knot_u,
	                               std::vector<double> &knot_v,
	                               std::vector<double> &knot_w) const;
	int nMeshRectangles()         const                { return meshrect_.size(); };

	// more get-methods
	std::vector<MeshRectangle*>::iterator  meshrectBegin()          { return meshrect_.begin(); };
	std::vector<MeshRectangle*>::iterator  meshrectEnd()            { return meshrect_.end(); };
	const std::vector<MeshRectangle*>& getAllMeshRectangles() const { return meshrect_ ;};
	MeshRectangle* getMeshRectangle(int i)                          { return meshrect_[i]; };
	const MeshRectangle* getMeshRectangle(int i) const              { return meshrect_[i]; };
	void getBezierElement(   int iEl, std::vector<double> &controlPoints) const;
	void getBezierExtraction(int iEl, std::vector<double> &extractMatrix) const;
	int getMaxContinuity(int i) const ;
	int getMinContinuity(int i) const ;
	//! \brief sets the maximum continuity in direction *dir to *cont. Note that this never _increases_ continuity, only decreases 
	void setMaxContinuity(int dir, int maxCont) ;

	// assorted specialized functions
	std::vector<Meshline*> getEdgeKnots(parameterEdge edge, bool normalized=false) const;
	void getDiagonalElements(      std::vector<int> &result) const ;
	void getDiagonalBasisfunctions(std::vector<int> &result) const ;
	void printElements(std::ostream &out) const;

	// fetch function spaces of different order/continuity
  LRSplineVolume*  getDerivedBasis(int raise_p1, int raise_p2, int raise_p3, size_t lower_k1, size_t lower_k2, size_t lower_k3, int dim=1) const;

	// interpolate and approximate functions
	/*
	void getLeastSquaresEdge(double (*f)(double, double),
	                         parameterEdge edge,
	                         std::vector<int> id,
	                         std::vector<double> val) const;
	*/

	//! \brief sets the controlpoints of *this to be a variational diminishing spline approximation of the input argument
	bool setControlPointsVDSA(const LRSplineVolume *lr);

	// input output methods
	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

	MeshRectangle* insert_line(MeshRectangle *newRect) ;

private:
	// make non-copyable
	LRSplineVolume(const LRSplineVolume&) = delete;
	LRSplineVolume& operator=(const LRSplineVolume&) = delete;

	// caching stuff
	mutable std::vector<std::vector<std::vector<int> > > elementCache_;
	mutable std::vector<double>            glob_knot_u_;
	mutable std::vector<double>            glob_knot_v_;
	mutable std::vector<double>            glob_knot_w_;
	mutable bool                           builtElementCache_;

	void createElementCache() const;

	std::vector<MeshRectangle*> meshrect_;

	void aPosterioriFixElements();
	void split(int constDir, Basisfunction *b, double new_knot, int multiplicity, HashSet<Basisfunction*> &newFunctions);
	void initMeta();
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3,
	          typename RandomIterator4>
	void initCore(int n1, int n2, int n3, int order_u, int order_v, int order_w, RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 knot_w, RandomIterator4 coef, int dim, bool rational=false) {
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

		int p1       = order_u;
		int p2       = order_v;
		int p3       = order_w;
		int unique_u = 0;
		int unique_v = 0;
		int unique_w = 0;
		std::vector<int> elm_u; // element index of the knot vector (duplicate knots is multiple index in knot_u, single index in elmRows)
		std::vector<int> elm_v;
		std::vector<int> elm_w;
		for(int i=0; i<n1+order_u; i++) {// const u, spanning v,w
			int mult = 1;
			elm_u.push_back(unique_u);
			while(i+1<n1+order_u && knot_u[i]==knot_u[i+1]) {
				i++;
				mult++;
			  elm_u.push_back(unique_u);
			}
			unique_u++;
			meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],  knot_w[0],
			                                      knot_u[i], knot_v[n2], knot_w[n3], mult));
		}
		for(int i=0; i<n2+order_v; i++) {// const v, spanning u,w
			int mult = 1;
			elm_v.push_back(unique_v);
			while(i+1<n2+order_v && knot_v[i]==knot_v[i+1]) {
				i++;
				mult++;
			  elm_v.push_back(unique_v);
			}
			unique_v++;
			meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[i], knot_w[0],
			                                      knot_u[n1], knot_v[i], knot_w[n3], mult));
		}
		for(int i=0; i<n3+order_w; i++) {
			int mult = 1;
			elm_w.push_back(unique_w);
			while(i+1<n3+order_w && knot_w[i]==knot_w[i+1]) {
				i++;
				mult++;
			  elm_w.push_back(unique_w);
			}
			unique_w++;
			meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[0],  knot_w[i],
			                                      knot_u[n1], knot_v[n2], knot_w[i], mult));
		}
		std::vector<std::vector<std::vector<Element*> > > elmRows(unique_w-1, std::vector<std::vector<Element*> >(unique_v-1));
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
					Element *elm = new Element(3, min, max);
					element_.push_back(elm);
					elmRows[k][j].push_back(elm);
				}
			}
		}

		for(int k=0; k<n3; k++)
			for(int j=0; j<n2; j++)
				for(int i=0; i<n1; i++) {
					RandomIterator1 ku = knot_u + i;
					RandomIterator2 kv = knot_v + j;
					RandomIterator3 kw = knot_w + k;
					RandomIterator4 c = coef + (k*n1*n2 + j*n1 + i)*(dim + rational);
					Basisfunction *b = new Basisfunction(ku, kv, kw, c , dim, order_u, order_v, order_w);
					basis_.insert(b);
					for(int k1=elm_w[k]; k1<elm_w[k+p3]; k1++)
						for(int k2=elm_v[j]; k2<elm_v[j+p2]; k2++)
							updateSupport(b, elmRows[k1][k2].begin() + elm_u[i], elmRows[k1][k2].begin() + elm_u[i+p1]);
		}
	}

};

} // end namespace LR

#endif

