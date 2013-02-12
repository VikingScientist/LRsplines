
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Basisfunction.h"

namespace LR {

Element::Element() {
	id_      = -1;
	overloadCount = 0;
}

Element::Element(int dim) {
	id_      = -1;
	overloadCount = 0;
	min.resize(dim);
	max.resize(dim);
}

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

void Element::removeSupportFunction(Basisfunction *f) {
	support_.erase(f);
}

void Element::addSupportFunction(Basisfunction *f) {
	support_.insert(f);
}

Element* Element::copy() {
	Element *returnvalue = new Element();
	
	returnvalue->id_          = this->id_;
	returnvalue->min          = this->min;    // it seems that the default vector operator= thing takes a deep copy 
	returnvalue->max          = this->max;
	returnvalue->support_ids_ = this->support_ids_;
	
	return returnvalue;
}


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
				someChange = true;
				break;
			}
		}
	}
	return newElement;
}

void Element::updateBasisPointers(std::vector<Basisfunction*> &basis) {
	for(uint i=0; i<support_ids_.size(); i++) {
		// add pointer from Element to Basisfunction
		support_.insert(basis[support_ids_[i]]);
		// add pointer from Basisfunction back to Element
		basis[support_ids_[i]]->addSupport(this);
	}
}

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

bool Element::isOverloaded()  const {
	int n = support_.size();
	if(n > 0) {
		int p1 = (*support_.begin())->getOrder(0);
		int p2 = (*support_.begin())->getOrder(1);
		if(n > p1*p2)
			return true;
	}
	return false;
}

/*
int Element::overloadedBasisCount() const {
	int ans = 0;
	for(uint i=0; i<support_.size(); i++)
		if(support_[i]->isOverloaded())
			ans++;
	return ans;
}
*/

} // end namespace LR
