
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

namespace LR {

Element::Element() {
	start_u_ =  0;
	start_v_ =  0;
	stop_u_  =  0;
	stop_v_  =  0;
	id_      = -1;
}

Element::Element(double start_u, double start_v, double stop_u, double stop_v) {

	start_u_ = start_u;
	start_v_ = start_v;
	stop_u_  = stop_u ;
	stop_v_  = stop_v ;
	id_      = -1;
}

void Element::removeSupportFunction(Basisfunction *f) {
	for(uint i=0; i<support_.size(); i++) {
		if(*f == *support_[i]) {
			support_[i] = support_.back();
			support_[support_.size()-1] = NULL;
			support_.pop_back();
			return;
		}
	}
}

void Element::addSupportFunction(Basisfunction *f) {
	for(uint i=0; i<support_.size(); i++) {
		if(f == support_[i]) {
			return;
		}
	}
	support_.push_back(f);
	f->addSupport(this);
}

Element* Element::split(bool split_u, double par_value) {
	Element *newElement = NULL;
	if(split_u) {
		if(par_value >= stop_u_ || par_value <= start_u_)
			return NULL;
		newElement = new Element(par_value, start_v_, stop_u_, stop_v_);
		stop_u_ = par_value;
	} else {
		if(par_value >= stop_v_ || par_value <= start_v_)
			return NULL;
		newElement = new Element(start_u_, par_value, stop_u_, stop_v_);
		stop_v_ = par_value;
	}
	for(uint i=0; i<support_.size(); i++) {
		if(support_[i]->overlaps(newElement))
			newElement->addSupportFunction(support_[i]);
		if(!support_[i]->overlaps(this)) {
			support_[i]->removeSupport(this);
			support_[i] = support_.back();
			support_.pop_back();
			i--;
		}
	}
	return newElement;
}

void Element::updateBasisPointers(std::vector<Basisfunction*> &basis) {
	for(uint i=0; i<support_ids_.size(); i++) {
		// add pointer from Element to Basisfunction
		support_.push_back(basis[support_ids_[i]]);
		// add pointer from Basisfunction back to Element
		support_.back()->addSupport(this);
	}
}

// convenience macro for reading formated input
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing element\n"; exit(326); } ws(is); }
void Element::read(std::istream &is) {
	char nextChar;
	is >> id_;
	ASSERT_NEXT_CHAR(':');
	ASSERT_NEXT_CHAR('(');
	is >> start_u_;
	ASSERT_NEXT_CHAR(',');
	is >> stop_u_;
	ASSERT_NEXT_CHAR(')');
	ASSERT_NEXT_CHAR('x');
	ASSERT_NEXT_CHAR('(');
	is >> start_v_;
	ASSERT_NEXT_CHAR(',');
	is >> stop_v_;
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
	os << id_ << ": (" << start_u_ << ", " << stop_u_ << ") x (" << start_v_ << ", " << stop_v_ << ")";
	os << "    {";
	for(uint i=0; i<support_.size()-1; i++) {
		os << support_[i]->getId() << ", " ;
	}
	os << support_.back()->getId() << "}";
}

} // end namespace LR
