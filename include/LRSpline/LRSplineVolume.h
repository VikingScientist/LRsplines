#ifndef LRSPLINEVOLUME_H
#define LRSPLINEVOLUME_H

#ifdef HAS_GOTOOLS
	#include <GoTools/utils/Point.h>
	#include <GoTools/trivariate/SplineVolume.h>
#endif
#ifdef HAS_BOOST
#include <boost/rational.hpp>
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
		initMeta();
		initCore(n1, n2, n3, order_u, order_v, order_w, knot_u, knot_v, knot_w, coef, dim, rational);
	}
	// LRSplineVolume(int n1, int n2, int n3, int order_u, int order_v, int order_w, double *knot_u, double *knot_v, double *knot_w, double *coef, int dim, bool rational=false);
	~LRSplineVolume();

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

	// linear independence methods
	bool isLinearIndepByOverloading(bool verbose) ;
	bool isLinearIndepByMappingMatrix(bool verbose) const { return true; };
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
	void getGlobalKnotUniqueVector(std::vector<double> &knot_u,
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
	void getDiagonalElements(      std::vector<int> &result) const ;
	void getDiagonalBasisfunctions(std::vector<Basisfunction*> &result) const ;
	void printElements(std::ostream &out) const;

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
	std::vector<MeshRectangle*> meshrect_;

	void split(int constDir, Basisfunction *b, double new_knot, int multiplicity, HashSet<Basisfunction*> &newFunctions);
	void initMeta();
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3,
	          typename RandomIterator4>
	void initCore(int n1, int n2, int n3, int order_u, int order_v, int order_w, RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 knot_w, RandomIterator4 coef, int dim, bool rational=false) {
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
	
		for(int k=0; k<n3; k++)
			for(int j=0; j<n2; j++)
				for (int i = 0; i < n1; i++) {
					RandomIterator1 ku = knot_u + i;
					RandomIterator2 kv = knot_v + j;
					RandomIterator3 kw = knot_w + k;
					RandomIterator4 c = coef + (k*n1*n2 + j*n1 + i)*(dim + rational);
					basis_.insert(new Basisfunction(ku, kv, kw, c , dim, order_u, order_v, order_w));
		}
		int unique_u=0;
		int unique_v=0;
		int unique_w=0;
		for(int i=0; i<n1+order_u; i++) {// const u, spanning v,w
			int mult = 1;
			while(i+1<n1+order_u && knot_u[i]==knot_u[i+1]) {
				i++;
				mult++;
			}
			unique_u++;
			meshrect_.push_back(new MeshRectangle(knot_u[i], knot_v[0],  knot_w[0],
			                                      knot_u[i], knot_v[n2], knot_w[n3], mult));
		}
		for(int i=0; i<n2+order_v; i++) {// const v, spanning u,w
			int mult = 1;
			while(i+1<n2+order_v && knot_v[i]==knot_v[i+1]) {
				i++;
				mult++;
			}
			unique_v++;
			meshrect_.push_back(new MeshRectangle(knot_u[0],  knot_v[i], knot_w[0],
			                                      knot_u[n1], knot_v[i], knot_w[n3], mult));
		}
		for(int i=0; i<n3+order_w; i++) {
			int mult = 1;
			while(i+1<n3+order_w && knot_w[i]==knot_w[i+1]) {
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
	

};

} // end namespace LR

#endif

