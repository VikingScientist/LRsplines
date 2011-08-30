
#include "Element.h"
#include "Basisfunction.h"

namespace LR {

Element::Element(double start_u, double start_v, double stop_u, double stop_v) {

	start_u_ = start_u;
	start_v_ = start_v;
	stop_u_  = stop_u ;
	stop_v_  = stop_v ;
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

void Element::read(std::istream &is) {
}

void Element::write(std::ostream &os) const {
	os << "(umin, umax) x (vmin vmax) = (" << start_u_ << ", " << stop_u_ << ") x (" << start_v_ << ", " << stop_v_ << ")\n";
	os << "support: {\n";
	for(uint i=0; i<support_.size(); i++) {
		os << "\t" << *support_[i] << std::endl;
	}
	os << "}";
}



} // end namespace LR
