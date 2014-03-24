#ifndef LRSPLINESURFACE_H
#define LRSPLINESURFACE_H

#include <vector>
#ifdef HAS_GOTOOLS
	#include <GoTools/utils/Point.h>
	#include <GoTools/geometry/SplineSurface.h>
#endif
#ifdef HAS_BOOST
#include <boost/rational.hpp>
#endif
#include "LRSpline.h"
#include "Basisfunction.h"
#include "Meshline.h"
#include "Element.h"
#include "HashSet.h"


namespace LR {


class LRSplineSurface : public LRSpline {

public:
	// constructors and destructors
	LRSplineSurface();
#ifdef HAS_GOTOOLS
	LRSplineSurface(Go::SplineSurface *surf);
#endif
	LRSplineSurface(int n1, int n2, int order_u, int order_v);
	LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v);
	LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational=false);
	~LRSplineSurface();
	LRSplineSurface* copy() const;
		
	// surface evaluation
#ifdef HAS_GOTOOLS
	virtual void point(Go::Point &pt, double u, double v, int iEl=-1) const;
	virtual void point(Go::Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const;
	virtual void point(std::vector<Go::Point> &pts, double upar, double vpar, int derivs, int iEl=-1) const;
	void computeBasis (double param_u, double param_v, Go::BasisPtsSf     & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf  & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf2 & result, int iEl=-1 ) const;
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

	// refinement functions
	void refineBasisFunction(int index);
	void refineBasisFunction(const std::vector<int> &indices);
	void refineElement(int index);
	void refineElement(const std::vector<int> &indices);
	void refineByDimensionIncrease(const std::vector<double> &error, double beta);

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

	// linear independence methods
	bool isLinearIndepByOverloading(bool verbose) ;
	bool isLinearIndepByFloatingPointMappingMatrix(bool verbose) const ;
#ifdef HAS_BOOST
	bool isLinearIndepByMappingMatrix(bool verbose) const ;
	void getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const ;
#endif

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

	// assorted specialized functions
	double makeIntegerKnots();
	void getSupportElements(       std::vector<int> &result, const std::vector<int> &basisfunctions) const ;
	void getDiagonalElements(      std::vector<int> &result) const ;
	void getDiagonalBasisfunctions(std::vector<int> &result) const ;
	void printElements(std::ostream &out) const;
	LRSplineSurface* getRaiseOrderSpace(int raiseOrderU, int raiseOrderV) const;
	std::vector<LRSplineSurface*> getDerivativeSpace() const ;
	bool setGlobalContinuity(int contU, int contV);
	bool decreaseContinuity( int du,    int dv);

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
	// initializeation methods (called from constructors)
	void initMeta();
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3>
	void initCore(int n1, int n2, int order_u, int order_v, RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 coef, int dim, bool rational=false) {
		int p1     = order_u;
		int p2     = order_v;
		order_.resize(2);
		start_.resize(2);
		end_.resize(2);
		rational_ = rational;
		dim_      = dim;
		order_[0]  = order_u;
		order_[1]  = order_v;
		start_[0]  = knot_u[0];
		start_[1]  = knot_v[0];
		end_[0]    = knot_u[n1+p1-1];
		end_[1]    = knot_v[n2+p2-1];
	
		for(int j=0; j<n2; j++)
			for(int i=0; i<n1; i++)
				basis_.insert(new Basisfunction(knot_u+i, knot_v+j, coef+(j*n1+i)*(dim+rational), dim, order_u, order_v));
		int unique_u=0;
		int unique_v=0;
		for(int i=0; i<n1+p1; i++) {// const u, spanning v
			int mult = 1;
			while(i<n1+p1 && knot_u[i]==knot_u[i+1]) {
				i++;
				mult++;
			}
			unique_u++;
			meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2+p2-1], mult) );
		}
		for(int i=0; i<n2+p2; i++) {// const v, spanning u
			int mult = 1;
			while(i<n2+p2 && knot_v[i]==knot_v[i+1]) {
				i++;
				mult++;
			}
			unique_v++;
			meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n1+p1-1], mult) );
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

		HashSet_iterator<Basisfunction*> it;
		for(it=basis_.begin(); it!=basis_.end(); ++it)
			updateSupport(*it);
	}

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

