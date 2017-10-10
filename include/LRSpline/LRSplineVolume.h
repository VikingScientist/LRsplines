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
	int getElementContaining(double u, double v, double w) const;
	// TODO: get rid of the iEl argument in evaluation signatures - it's too easy to mess it up (especially with derivatives at multiple-knot boundaries).
	//       Try and sort the Elements after all refinements and binary search for the containing point in logarithmic time

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

	// input output methods
	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

	MeshRectangle* insert_line(MeshRectangle *newRect) ;

private:
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
		int p1     = order_u;
		int p2     = order_v;
		int p3     = order_w;
		// sanity check input
		if(n1 < p1 ||
		   n2 < p2 ||
		   n3 < p3) {
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
		order_[0]  = p1;
		order_[1]  = p2;
		order_[2]  = p3;
		start_[0]  = knot_u[p1-1];
		start_[1]  = knot_v[p2-1];
		start_[2]  = knot_w[p3-1];
		end_[0]    = knot_u[n1];
		end_[1]    = knot_v[n2];
		end_[2]    = knot_w[n3];

		std::vector<double> elm_u; // list of inner knots used to create elements (ingnoring ghost domain)
		std::vector<double> elm_v;
		std::vector<double> elm_w;
		std::vector<int> elm_u_i;
		std::vector<int> elm_v_i;
		std::vector<int> elm_w_i;
		for(int i=0; i<n1+p1; i++) {// const u, spanning v,w
			int mult = 1;
			elm_u_i.push_back(elm_u.size());
			while(i+1<n1+p1 && knot_u[i]==knot_u[i+1]) {
				i++;
				mult++;
			  elm_u_i.push_back(elm_u.size());
			}
			if(p1-1<=i && i-mult<n1)
				elm_u.push_back(knot_u[i]);
			if(i-mult >= n1)
				elm_u_i.back()--;
			meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],       knot_w[0],
			                                      knot_u[i], knot_v[n2+p2-1], knot_w[n3+p3-1], mult));
		}
		for(int i=0; i<n2+p2; i++) {// const v, spanning u,w
			int mult = 1;
			elm_v_i.push_back(elm_v.size());
			while(i+1<n2+p2 && knot_v[i]==knot_v[i+1]) {
				i++;
				mult++;
			  elm_v_i.push_back(elm_v.size());
			}
			if(p2-1<=i && i-mult<n2)
				elm_v.push_back(knot_v[i]);
			if(i-mult >= n2)
				elm_v_i.back()--;
			meshrect_.push_back(new MeshRectangle(knot_u[0],       knot_v[i], knot_w[0],
			                                      knot_u[n1+p1-1], knot_v[i], knot_w[n3+p3-1], mult));
		}
		for(int i=0; i<n3+p3; i++) {
			int mult = 1;
			elm_w_i.push_back(elm_w.size());
			while(i+1<n3+p3 && knot_w[i]==knot_w[i+1]) {
				i++;
				mult++;
			  elm_w_i.push_back(elm_w.size());
			}
			if(p3-1<=i && i-mult<n3)
				elm_w.push_back(knot_w[i]);
			if(i-mult >= n3)
				elm_w_i.back()--;
			meshrect_.push_back(new MeshRectangle(knot_u[0],       knot_v[0],       knot_w[i],
			                                      knot_u[n1+p1-1], knot_v[n2+p2-1], knot_w[i], mult));
		}
		std::vector<std::vector<std::vector<Element*> > > elmRows(elm_w.size()-1, std::vector<std::vector<Element*> >(elm_v.size()-1));
		for(int k=0; k<elm_w.size()-1; k++) {
			for(int j=0; j<elm_v.size()-1; j++) {
				for(int i=0; i<elm_u.size()-1; i++) {
					double umin = elm_u[ i ];
					double umax = elm_u[i+1];
					double vmin = elm_v[ i ];
					double vmax = elm_v[i+1];
					double wmin = elm_w[ i ];
					double wmax = elm_w[i+1];
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
					Basisfunction *b = new Basisfunction(ku, kv, kw, c , dim, p1, p2, p3);
					basis_.insert(b);
					for(int k1=elm_w_i[k]; k1<elm_w_i[k+p3]; k1++)
						for(int k2=elm_v_i[j]; k2<elm_v_i[j+p2]; k2++)
							updateSupport(b, elmRows[k1][k2].begin() + elm_u_i[i], elmRows[k1][k2].begin() + elm_u_i[i+p1]);
		}
	}

};

} // end namespace LR

#endif

