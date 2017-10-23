
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Basisfunction.h"
#include <stdlib.h>

typedef unsigned int uint;

namespace LR {

/************************************************************************************************************************//**
 * \brief Default constructor
 ***************************************************************************************************************************/
Element::Element() {
	id_      = -1;
	overloadCount = 0;
}

/************************************************************************************************************************//**
 * \brief Constructor
 * \param dim The dimension of the element (2 for SplineSurfaces, 3 for SplineVolumes)
 ***************************************************************************************************************************/
Element::Element(int dim) {
	id_           = -1;
	overloadCount = 0;
	min.resize(dim);
	max.resize(dim);
}

/************************************************************************************************************************//**
 * \brief Surface element constructor
 * \param start_u Lower left u-coordinate
 * \param start_v Lower left v-coordinate
 * \param stop_u  Upper right u-coordinate
 * \param stop_v  Upper right v-coordinate
 ***************************************************************************************************************************/
Element::Element(double start_u, double start_v, double stop_u, double stop_v) {
	min.resize(2);
	max.resize(2);

	min[0] = start_u;
	min[1] = start_v;
	max[0] = stop_u ;
	max[1] = stop_v ;
	id_    = -1;
	overloadCount = 0;
}

/************************************************************************************************************************//**
 * \brief Removes one basisfunction from the list of supported Basisfunctions
 * \param f The pointer to the Basisfunction to remove
 ***************************************************************************************************************************/
void Element::removeSupportFunction(Basisfunction *f) {
	support_.erase(f);
}

/************************************************************************************************************************//**
 * \brief Adds one basisfunction from the list of supported Basisfunctions
 * \param f The pointer to the Basisfunction to add
 ***************************************************************************************************************************/
void Element::addSupportFunction(Basisfunction *f) {
	support_.insert(f);
}

/************************************************************************************************************************//**
 * \brief Makes a deep copy of the Element and returns a pointer to this (caller responsible for freeing memory). Note: Does not copy support Basisfunction
 ***************************************************************************************************************************/
Element* Element::copy() {
	Element *returnvalue = new Element();

	returnvalue->id_          = this->id_;
	returnvalue->min          = this->min;    // it seems that the default vector operator= thing takes a deep copy
	returnvalue->max          = this->max;

	for(Basisfunction* b : support_)
		returnvalue->support_ids_.push_back(b->getId());

	return returnvalue;
}


/************************************************************************************************************************//**
 * \brief Splits an element into two new ones by reducing the size of this element and returning a new element.
 * \param splitDim The constant parameter direction to split the element
 * \param par_value The parameter value to do the splitting
 * \returns The new element resulting from the splitting
 ***************************************************************************************************************************/
Element* Element::split(int splitDim, double par_value) {
	Element *newElement = NULL;
	if(par_value >= max[splitDim] || par_value <= min[splitDim])
		return NULL;

	std::vector<double> newMin(min.begin(), min.end());
	std::vector<double> newMax(max.begin(), max.end());

	newMin[splitDim] = par_value; // new element should start at par_value
	max[splitDim]    = par_value; // old element should stop  at par_value

	newElement = new Element(min.size(), newMin.begin(), newMax.begin());

	for(Basisfunction *b : support_)
		if(b->addSupport(newElement)) // tests for overlapping as well
			newElement->addSupportFunction(b);

	bool someChange = true;
	while(someChange) {
		someChange = false;
		for(Basisfunction *b : support_) {
			if(!b->overlaps(this)) {
				support_.erase(b);
				b->removeSupport(this);
				someChange = true;
				break;
			}
		}
	}
	return newElement;
}

/************************************************************************************************************************//**
 * \brief returns the parametric midpoint of the element, i.e. [ (umin()+umax())/2, (vmin()+vmax())/2] for 2D elements
 ***************************************************************************************************************************/
std::vector<double> Element::midpoint() const {
	std::vector<double> result;
	for(uint i=0; i<min.size(); i++)
		result.push_back((min[i]+max[i])/2.0);
	return result;
}

/************************************************************************************************************************//**
 * \brief Only used at the end of the read function. Updates the element-to-basisfunction and back again pointers.
 * \param basis The flat vector list of basisfunctions
 ***************************************************************************************************************************/
void Element::updateBasisPointers(std::vector<Basisfunction*> &basis) {
	for(uint i=0; i<support_ids_.size(); i++) {
		// add pointer from Element to Basisfunction
		support_.insert(basis[support_ids_[i]]);
		// add pointer from Basisfunction back to Element
		basis[support_ids_[i]]->addSupport(this);
	}
}

/************************************************************************************************************************//**
 * \brief Reads formatted input from input stream
 * \param is The input stream to read from
 ***************************************************************************************************************************/
// convenience macro for reading formated input
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing element\n"; exit(326); } ws(is); }
void Element::read(std::istream &is) {
	char nextChar;
	int dim;
	is >> id_;
	ASSERT_NEXT_CHAR('[');
	is >> dim;
	ASSERT_NEXT_CHAR(']');
	ASSERT_NEXT_CHAR(':');
	min.resize(dim);
	max.resize(dim);

	ASSERT_NEXT_CHAR('(');
	is >> min[0];
	for(int i=1; i<dim; i++) {
		ASSERT_NEXT_CHAR(',');
		is >> min[i];
	}
	ASSERT_NEXT_CHAR(')');
	ASSERT_NEXT_CHAR('x');
	ASSERT_NEXT_CHAR('(');
	is >> max[0];
	for(int i=1; i<dim; i++) {
		ASSERT_NEXT_CHAR(',');
		is >> max[i];
	}
	ASSERT_NEXT_CHAR(')');
	ASSERT_NEXT_CHAR('{');

	// read id's of all supported basis functions
	int basis_id;
	is >> basis_id;
	support_ids_.push_back(basis_id);
	ws(is);
	nextChar = is.peek();
	while(nextChar == ',') {
		is.get(); ws(is);
		is >> basis_id;
		support_ids_.push_back(basis_id);
		nextChar = is.peek();
	}
	ASSERT_NEXT_CHAR('}');
}
#undef ASSERT_NEXT_CHAR

/************************************************************************************************************************//**
 * \brief Writes the element primitive to output stream
 * \param os The output stream to write to
 ***************************************************************************************************************************/
void Element::write(std::ostream &os) const {
	os << id_ << " [" << min.size() << "] : ";
	os << "(" << min[0];
	for(uint i=1; i<min.size(); i++)
		os << ", " << min[i] ;
	os << ") x (" << max[0];
	for(uint i=1; i<max.size(); i++)
		os << ", " << max[i] ;
	os << ")";
	os << "    {";
	bool isFirst = true;
	for(Basisfunction *b : support_) {
		if(!isFirst) os << ", ";
		os << b->getId();
		isFirst = false;
	}
	os << "}";
}

/************************************************************************************************************************//**
 * \brief Checks if an element is overloaded by having too many Basisfunctions with support on this element
 * \returns True if more than (p+1)*(q+1) Basisfunctions have support on this element, where p and q are the polynomial degrees.
            In case of trivariate LRspline Volumes (p+1)*(q+1)*(r+1) is used
 ***************************************************************************************************************************/
bool Element::isOverloaded()  const {
	int n = support_.size();
	Basisfunction *b = *support_.begin();
	if(n > 0) {
		if(b->nVariate() == 2) { // surfaces
			int p1 = (*support_.begin())->getOrder(0);
			int p2 = (*support_.begin())->getOrder(1);
			if(n > p1*p2)
				return true;
		} else if(b->nVariate() == 3) { // volumes
			int p1 = (*support_.begin())->getOrder(0);
			int p2 = (*support_.begin())->getOrder(1);
			int p3 = (*support_.begin())->getOrder(2);
			if(n > p1*p2*p3)
				return true;
		}
	}
	return false;
}

} // end namespace LR
