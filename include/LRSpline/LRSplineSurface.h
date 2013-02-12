#ifndef LR_SPLINE_H
#define LR_SPLINE_H

#include <vector>
#include <set>
#include <GoTools/utils/Point.h>
#include <GoTools/geometry/Streamable.h>
#include <boost/rational.hpp>
#include <GoTools/geometry/SplineSurface.h>
#include "Basisfunction.h"
#include "HashSet.h"

enum refinementStrategy {
	LR_SAFE = 0,
	LR_MINSPAN = 1,
	LR_ISOTROPIC_EL = 2,
	LR_ISOTROPIC_FUNC = 3,
};


namespace LR {

class Basisfunction;
class Meshline;
class Element;

class LRSplineSurface : public Go::Streamable {

public:
	LRSplineSurface();
	LRSplineSurface(Go::SplineSurface *surf);
	LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational=false);
	~LRSplineSurface();

	LRSplineSurface &operator=(LRSplineSurface& copythis);
	LRSplineSurface* copy();
	// surface evaluation
	virtual void point(Go::Point &pt, double u, double v, int iEl=-1) const;
	virtual void point(Go::Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const;
	virtual void point(std::vector<Go::Point> &pts, double upar, double vpar, int derivs, int iEl=-1) const;
	void computeBasis (double param_u, double param_v, Go::BasisPtsSf     & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf  & result, int iEl=-1 ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf2 & result, int iEl=-1 ) const;
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

	// set refinement state parameters
	void setRefStrat(enum refinementStrategy strat) { refStrat_        = strat;    };
	void setRefSymmetry(int symmetry)               { this->symmetry_  = symmetry; };
	void setRefMultiplicity(int mult)               { refKnotlineMult_ = mult;     };
	void setMaxTjoints(int n)                       { maxTjoints_      = n;        };
	void setCloseGaps(bool doClose)                 { doCloseGaps_     = doClose;  };
	void setMaxAspectRatio(double r, bool aposterioriFix=true) {
		maxAspectRatio_ = r;
		doAspectRatioFix_ = aposterioriFix;
	}

	// (private) refinement functions
	Meshline* insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity=1);
	Meshline* insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity=1);
	void getFullspanLines(  int iEl,    std::vector<Meshline*>& lines);
	void getMinspanLines(   int iEl,    std::vector<Meshline*>& lines);
	void getStructMeshLines(int iBasis, std::vector<Meshline*>& lines);
	void aPosterioriFixes() ;
	void closeGaps(            std::vector<Meshline*>* newLines=NULL);
	void enforceMaxTjoints(    std::vector<Meshline*>* newLines=NULL);
	void enforceMaxAspectRatio(std::vector<Meshline*>* newLines=NULL);

	// linear independence methods
	bool isLinearIndepByOverloading(bool verbose) ;
	bool isLinearIndepByMappingMatrix(bool verbose) const ;
	bool isLinearIndepByFloatingPointMappingMatrix(bool verbose) const ;
	void getNullSpace(std::vector<std::vector<boost::rational<long long> > >& nullspace) const ;

	void updateSupport(Basisfunction *f) ;
	void updateSupport(Basisfunction *f,
	                   std::vector<Element*>::iterator start,
	                   std::vector<Element*>::iterator end ) ;
	
	void generateIDs() const;

	// common get methods
	void getGlobalKnotVector      (std::vector<double> &knot_u, std::vector<double> &knot_v) const;
	void getGlobalUniqueKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const;
	//void change_coefs(std::vector<double> coefs);
	virtual double startparam_u() const                { return start_u_; };
	virtual double startparam_v() const                { return start_v_; };
	virtual double endparam_u()   const                { return end_u_; };
	virtual double endparam_v()   const                { return end_v_; };
	virtual int dimension()       const                { return dim_; };
	int order_u()                 const                { return order_u_; };
	int order_v()                 const                { return order_v_; };
	bool rational()               const                { return rational_; };
	int nBasisFunctions()         const                { return basis_.size(); };
	int nElements()               const                { return element_.size(); };
	int nMeshlines()              const                { return meshline_.size(); };

	// more get-methods
	std::vector<Meshline*>::iterator meshlineBegin()               { return meshline_.begin(); };
	std::vector<Meshline*>::iterator meshlineEnd()                 { return meshline_.end(); };
	std::vector<Element*>::iterator elementBegin()                 { return element_.begin(); };
	std::vector<Element*>::iterator elementEnd()                   { return element_.end(); };
	HashSet_iterator<Basisfunction*>       basisBegin()            { return basis_.begin(); };
	HashSet_iterator<Basisfunction*>       basisEnd()              { return basis_.end(); };
	HashSet_const_iterator<Basisfunction*> basisBegin() const      { return basis_.begin(); };
	HashSet_const_iterator<Basisfunction*> basisEnd()   const      { return basis_.end(); };
	HashSet<Basisfunction*> getAllBasisfunctions()                 { return basis_ ;};
	Element* getElement(int i)                                     { return element_[i]; };
	Basisfunction* getBasisfunction(int iBasis) {
		HashSet_iterator<Basisfunction*> it = basis_.begin();
		for(int i=0; i<iBasis; i++) it++;
		return *it;
	}
	void getEdgeFunctions(std::vector<Basisfunction*> &edgeFunctions, parameterEdge edge, int depth=1) const;
	

	// assorted specialized functions
	void set_dim(int dimvalue)                         {dim_ = dimvalue;};
	void rebuildDimension(int dimvalue) ;
	double makeIntegerKnots();
	void getSupportElements(       std::vector<int> &result, const std::vector<int> &basisfunctions) const ;
	void getDiagonalElements(      std::vector<int> &result) const ;
	void getDiagonalBasisfunctions(std::vector<int> &result) const ;
	void printElements(std::ostream &out) const;

	// interpolate and approximate functions
	void getLeastSquaresEdge(double (*f)(double, double),
	                         parameterEdge edge,
	                         std::vector<int> id,
	                         std::vector<double> val) const;

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
	int split(bool insert_in_u, Basisfunction *b, double new_knot, int multiplicity=1);
	Meshline* insert_line(bool const_u, double const_par, double start, double stop, int multiplicity);
	
	bool rational_;
	HashSet<Basisfunction*> basis_;
	std::vector<Meshline*> meshline_;
	std::vector<Element*> element_;
	int dim_;
	int order_u_;
	int order_v_;
	double start_u_;
	double start_v_;
	double end_u_;
	double end_v_;

	// refinement parameters
	enum refinementStrategy refStrat_;
	int                     refKnotlineMult_;
	int                     symmetry_;
	int                     maxTjoints_;
	bool                    doCloseGaps_;
	bool                    doAspectRatioFix_;
	double                  maxAspectRatio_;

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

