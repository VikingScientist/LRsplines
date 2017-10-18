#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

#ifdef HAS_GOTOOLS
	#include <GoTools/utils/Point.h>
#endif
#include <vector>
#include "HashSet.h"
#include "Streamable.h"
#include "LRSpline.h"

namespace LR {

class Element;

/************************************************************************************************************************//**
 * \brief Basisfunction class to store the individual B-splines which make up the LR B-spline space
 * \details Stores the local knot vectors corresponding in each parametric direction (two for bivariate surfaces, three for
 *          trivariate volumes) as well as the control point and scaling weight. Used for evaluation of the B-splines and
 *          all their derivatives. The class does also have pointers back to the elements which they have support on
 ***************************************************************************************************************************/
class Basisfunction : public Streamable {
public:
	Basisfunction(int dim, int order_u, int order_v);
	/************************************************************************************************************************//**
	 * \brief Constructor for arbitray high parametric dimension
	 * \param physDim The dimension in the physical space, i.e. the number of components of the controlpoints
	 * \param parDim The dimension in the parametric space, i.e. the number of local knot vectors
	 * \param order List of polynomial orders (degree + 1) in each parametric direction
	 ***************************************************************************************************************************/
	template <typename RandomIterator>
	Basisfunction(int physDim, int parDim, RandomIterator order) {
		weight_       = 1;
		id_           = -1;
		hashCode_     = 0;
		knots_.resize(parDim);
		for(int i=0; i<parDim; i++)
			knots_[i].resize(order[i]+1);
		controlpoint_.resize(physDim);
	}

	/************************************************************************************************************************//**
	 * \brief Constructor for bivariate Basisfunctions
	 * \param knot_u Knot vector in first parametric direction
	 * \param knot_v Knot vector in second parametric direction
	 * \param controlpoint The control point associated with the B-spline
	 * \param dim Physical dimension, i.e. the number of components of the controlpoint
	 * \param order_u Polynomial order (degree + 1) in first parametric direction
	 * \param order_v Polynomial order (degree + 1) in second parametric direction
	 * \param weight Scaling weight for partition of unity (note: not NURBS rational weight)
	 ***************************************************************************************************************************/
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3>
	Basisfunction(RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 controlpoint, int dim, int order_u, int order_v, double weight=1.0) {
		weight_       = weight ;
		id_           = -1;
		hashCode_     = 0;
		knots_.resize(2);
		knots_[0].resize(order_u+1);
		knots_[1].resize(order_v+1);
		controlpoint_.resize(dim);

		std::copy(knot_u,       knot_u       + order_u+1,   knots_[0].begin());
		std::copy(knot_v,       knot_v       + order_v+1,   knots_[1].begin());
		std::copy(controlpoint, controlpoint + dim,         controlpoint_.begin());
	}

	/************************************************************************************************************************//**
	 * \brief Constructor for trivariate Basisfunctions
	 * \param knot_u Knot vector in first parametric direction
	 * \param knot_v Knot vector in second parametric direction
	 * \param knot_w Knot vector in second parametric direction
	 * \param controlpoint The control point associated with the B-spline
	 * \param dim Physical dimension, i.e. the number of components of the controlpoint
	 * \param order_u Polynomial order (degree + 1) in first parametric direction
	 * \param order_v Polynomial order (degree + 1) in second parametric direction
	 * \param order_w Polynomial order (degree + 1) in third parametric direction
	 * \param weight Scaling weight for partition of unity (note: not NURBS rational weight)
	 ***************************************************************************************************************************/
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3,
	          typename RandomIterator4>
	Basisfunction(RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 knot_w, RandomIterator4 controlpoint, int dim, int order_u, int order_v, int order_w, double weight=1.0) {
		weight_       = weight ;
		id_           = -1;
		hashCode_     = 0;
		knots_.resize(3);
		knots_[0].resize(order_u+1);
		knots_[1].resize(order_v+1);
		knots_[2].resize(order_w+1);
		controlpoint_.resize(dim);

		std::copy(knot_u,       knot_u       + order_u+1,   knots_[0].begin());
		std::copy(knot_v,       knot_v       + order_v+1,   knots_[1].begin());
		std::copy(knot_w,       knot_w       + order_w+1,   knots_[2].begin());
		std::copy(controlpoint, controlpoint + dim,         controlpoint_.begin());
	}
	virtual ~Basisfunction();
	Basisfunction* copy() const;

	//evaluation functions
	double evaluate(double u, double v, bool u_from_right=true, bool v_from_right=true) const;
	double evaluate(double u, double v, double w, bool u_from_right=true, bool v_from_right=true, bool w_from_right=true) const;
	void   evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right=true, bool v_from_right=true) const;
	void   evaluate(std::vector<double> &results, double u, double v, double w, int derivs, bool u_from_right=true, bool v_from_right=true, bool w_from_right=true) const;
	void   evaluate(std::vector<double> &results, const std::vector<double> &parPt, int derivs, const std::vector<bool> &from_right) const;

	// Basisfunction -> Element interatcion (support)
	bool                            overlaps(Element *el) const ;
	bool                            addSupport(Element *el)     ;
	bool                            removeSupport(Element *el)  ;
	std::vector<Element*>::iterator supportedElementBegin(){ return support_.begin(); };
	std::vector<Element*>::iterator supportedElementEnd()  { return support_.end();   };
	const std::vector<Element*>     support() const        { return support_;         };
	std::vector<Element*>           getExtendedSupport()        ;
	std::vector<Element*>           getMinimalExtendedSupport() ;
	HashSet<Basisfunction*>         getOverlappingFunctions() const ;
	bool                            isOverloaded()     const    ;
	int                             getOverloadCount() const    ;

	// get/set methods
	void setId(int id)  { this->id_ = id; };
	void setDimension(int dim)  ;

#ifdef HAS_GOTOOLS
	void      getControlPoint(Go::Point &pt) const;
	Go::Point getGrevilleParameter()         const;
#endif
	void   getControlPoint(std::vector<double> &pt)      const;
	void   getGrevilleParameter(std::vector<double> &pt) const;
	int    getId()                           const { return id_; };
	int    nSupportedElements()              const { return support_.size(); };
	int    nVariate()                        const { return knots_.size(); };
	int    dim()                             const { return controlpoint_.size(); };
	double getParmin(int i)                  const { return knots_[i][0];     };
	double getParmax(int i)                  const { return knots_[i].back(); };
	int    getOrder( int i)                  const { return knots_[i].size()-1;   };
	std::vector<double>& getknots(int i)           { return knots_[i]; };
	const std::vector<double>& getknots(int i) const { return knots_[i]; };
	std::vector<double>::iterator cp()             { return controlpoint_.begin(); };
	std::vector<double>::const_iterator cp() const { return controlpoint_.begin(); };
	double cp(int i)                         const { return controlpoint_[i]; };
	double w()                               const { return weight_; };

	long hashCode() const ;

	bool contains(const Basisfunction &other) const ;

	// operator overloading
	bool equals(const Basisfunction &other) const ;
	bool operator==(const Basisfunction &other) const;
	void operator+=(const Basisfunction &other) ;
	std::vector<double>&       operator[](int i)       { return knots_[i]; } ;
	const std::vector<double>& operator[](int i) const { return knots_[i]; } ;

	// IO-functions
	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

	void flip(int dir1=0, int dir2=1);
	void reverse(int pardir, double parmin=0.0, double parmax=1.0);
	void normalize(int pardir, double parmin, double parmax);

private:

	int                               id_;
	double                            weight_;
	mutable long                      hashCode_;
	std::vector<double>               controlpoint_;
	std::vector<std::vector<double> > knots_;
	std::vector<Element*>             support_;

};

} // end namespace LR

#endif

