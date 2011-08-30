
#include "Meshline.h"
#include "Element.h"
#include "Basisfunction.h"

namespace LR {

Meshline::Meshline(bool span_u_line, double const_par, double start, double stop, int multiplicity) {
	span_u_line_  =  span_u_line   ;
	const_par_    =  const_par     ;
	start_        =  start         ;
	stop_         =  stop          ;
	multiplicity_ =  multiplicity  ;
}

bool Meshline::containedIn(Basisfunction *basis) const {
	if(span_u_line_) {
		for(int i=0; i<=basis->order_v_; i++)
			if(basis->knot_v_[i] == const_par_)
				return true;
	} else { // span-v_line
		for(int i=0; i<=basis->order_u_; i++)
			if(basis->knot_u_[i] == const_par_)
				return true;
	}
	return false;
}

bool Meshline::splits(Element *el) const {
	if(span_u_line_) {
		if( el->vmin() < const_par_ && const_par_ < el->vmax() &&
		    start_ <= el->umin()    && el->umax() <= stop_)
			return  true;
	} else { // span-v line
		if( el->umin() < const_par_ && const_par_ < el->umax() &&
		    start_ <= el->vmin()    && el->vmax() <= stop_)
			return  true;
	}
	return false;
}

bool Meshline::splits(Basisfunction *basis) const {
	if(span_u_line_) {
		if( basis->knot_v_[0] < const_par_ && const_par_ < basis->knot_v_[basis->order_v_] &&
		    start_ <= basis->knot_u_[0]  && basis->knot_u_[basis->order_u_] <= stop_)
			return  true;
	} else { // span-v line
		if( basis->knot_u_[0] < const_par_ && const_par_ < basis->knot_u_[basis->order_u_] &&
		    start_ <= basis->knot_v_[0]  && basis->knot_v_[basis->order_v_] <= stop_)
			return  true;
	}
	return false;
}

bool Meshline::is_spanning_u() const {
	return span_u_line_;
}

void Meshline::read(std::istream &is) {
}

void Meshline::write(std::ostream &os) const {
	if(span_u_line_) 
		os <<  "[" << start_ << ", " << stop_ << "] x " << const_par_ << " (" << multiplicity_ << ")";
	else // span-v line
		os << const_par_ << " x [" << start_ << ", " << stop_ << "] (" << multiplicity_ << ")";
}


} // end namespace LR


