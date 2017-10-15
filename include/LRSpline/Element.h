#ifndef ELEMENT_H
#define ELEMENT_H

#include "Streamable.h"
#include <vector>
#include "HashSet.h"

namespace LR {

class Basisfunction;
class Meshline;

/************************************************************************************************************************//**
 * \brief Element class to partition the parametric space into subrectangles where all Basisfunctions are infitely differentiable
 * \details Stores the parametric bounding box of an element as well as a pointer to all the Basisfunctions which are active
 *          on this element. It is noteworthy to state that all computations on the Element class take place in the parametric
 *          space rather than in the physical (geometry) space. This class is shared by both LRSplineVolume and LRSplineSurface
 ***************************************************************************************************************************/
class Element : public Streamable {

public:
	Element();
	explicit Element(int dim);
	Element(double start_u, double start_v, double stop_u, double stop_v);
	/************************************************************************************************************************//**
	 * \brief General constructor
	 * \param dim The dimension of the element (i.e. 2 for surfaces, 3 for volumes)
	 * \param lowerLeft Template iterator to the lower left corner (for instance std::vector<double>::begin() or double*)
	 * \param upperRight Template iterator to the upper right corner
	 ***************************************************************************************************************************/
	template <typename RandomIterator1,
	          typename RandomIterator2>
	Element(int dim, RandomIterator1 lowerLeft, RandomIterator2 upperRight) {
		min.resize(dim);
		max.resize(dim);
		std::copy(lowerLeft,  lowerLeft  + dim, min.begin());
		std::copy(upperRight, upperRight + dim, max.begin());
		id_           = -1;
		overloadCount = 0;
	}
	Element(std::vector<double> &lowerLeft, std::vector<double> &upperRight);
	void removeSupportFunction(Basisfunction *f);
	void addSupportFunction(Basisfunction *f);
	Element *split(int splitDim, double par_value);
	Element* copy();
	virtual ~Element() {}
	// get/set methods
	//! \brief Get coordinate i of the lower left corner of the element
	double getParmin(int i) const { return min[i]; };
	//! \brief Get coordinate i of the upper right corner of the element
	double getParmax(int i) const { return max[i]; };
	double umin()           const { return min[0]; };
	double vmin()           const { return min[1]; };
	double wmin()           const { return min[2]; };
	double umax()           const { return max[0]; };
	double vmax()           const { return max[1]; };
	double wmax()           const { return max[2]; };
	std::vector<double> midpoint() const;
	//! \brief Returns the parametric area of the element
	double area()           const { return (max[1]-min[1])*(max[0]-min[0]);                  };
	//! \brief Returns the parametric volume of the element
	double volume()         const { return (max[2]-min[2])*(max[1]-min[1])*(max[0]-min[0]);  };
	HashSet_iterator<Basisfunction*> supportBegin()                 { return support_.begin(); };
	HashSet_iterator<Basisfunction*> supportEnd()                   { return support_.end();   };
	HashSet_const_iterator<Basisfunction*> constSupportBegin()const { return support_.begin(); };
	HashSet_const_iterator<Basisfunction*> constSupportEnd()  const { return support_.end();   };
	const HashSet<Basisfunction*>& support()                  const { return support_;   };
	//! \brief Returns the number of Basisfunctions with support on this element
	int nBasisFunctions() const           { return support_.size(); };
	//! \brief Sets a general purpose indexing id to this element
	void setId(int id)                    { this->id_ = id; };
	//! \brief Gets the id set by the Element::setId function
	int  getId() const                    { return id_; };
	//! \brief Gets the dimension of the element (2 for surfaces, 3 for volumes)
	int  getDim() const                   { return min.size(); };
	void setUmin(double u)                { min[0] = u; };
	void setVmin(double v)                { min[1] = v; };
	void setUmax(double u)                { max[0] = u; };
	void setVmax(double v)                { max[1] = v; };

	bool isOverloaded() const;
	void resetOverloadCount()    { overloadCount = 0;      }
	int incrementOverloadCount() { return overloadCount++; }
	int getOverloadCount() const { return overloadCount;   }

	void updateBasisPointers(std::vector<Basisfunction*> &basis) ;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

private:
	std::vector<double> min;  // lower left corner in typical 2 or 3 dimensions
	std::vector<double> max;  // upper right corner
	int id_;

	HashSet<Basisfunction*> support_;
	std::vector<int> support_ids_; // temporary storage for the read() method only

	int overloadCount ;

};

} // end namespace LR

#endif

