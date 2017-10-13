
#include "LRSpline/Meshline.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"
#include <cmath>
#include <cstdlib>

namespace LR {

#define DOUBLE_TOL 1e-14

Meshline::Meshline() {
	span_u_line_  = false;
	const_par_    = 0;
	start_        = 0;
	stop_         = 0;
	multiplicity_ = 0;
	type_         = INITIAL;
}

Meshline::Meshline(bool span_u_line, double const_par, double start, double stop, int multiplicity) {
	span_u_line_  =  span_u_line   ;
	const_par_    =  const_par     ;
	start_        =  start         ;
	stop_         =  stop          ;
	multiplicity_ =  multiplicity  ;
	type_         =  INITIAL       ;
}

Meshline::~Meshline() {
}


Meshline* Meshline::copy() {
	 Meshline *returnvalue     = new Meshline();

	 returnvalue->span_u_line_ = this->span_u_line_;
	 returnvalue->const_par_   = this->const_par_;
	 returnvalue->start_       = this->start_;
	 returnvalue->stop_        = this->stop_;
	 returnvalue->multiplicity_= this->multiplicity_;
	 returnvalue->type_        = this->type_;
	 return returnvalue;
}

int Meshline::nKnotsIn(Basisfunction *basis) const {
	int hits = 0;
	if(span_u_line_) {
		for(int i=0; i<=basis->getOrder(1); i++)
			if( fabs((*basis)[1][i] - const_par_) < DOUBLE_TOL )
				hits++;
	} else { // span-v_line
		for(int i=0; i<=basis->getOrder(0); i++)
			if( fabs((*basis)[0][i] - const_par_) < DOUBLE_TOL)
				hits++;
	}
	return hits;
}

bool Meshline::touches(Element *el) const {
	if(span_u_line_) {
		if( el->vmin() < const_par_ && const_par_ < el->vmax() &&
		   (start_ == el->umax()    || el->umin() == stop_))
			return  true;
	} else { // span-v line
		if( el->umin() < const_par_ && const_par_ < el->umax() &&
		   (start_ == el->vmax()    || el->vmin() == stop_))
			return  true;
	}
	return false;
}

bool Meshline::splits(Element *el) const {
	if(span_u_line_) {
		if( el->vmin() < const_par_ && const_par_ < el->vmax() &&
		    start_-DOUBLE_TOL <= el->umin()    && el->umax()-DOUBLE_TOL <= stop_)
			return  true;
	} else { // span-v line
		if( el->umin() < const_par_ && const_par_ < el->umax() &&
		    start_-DOUBLE_TOL <= el->vmin()    && el->vmax()-DOUBLE_TOL <= stop_)
			return  true;
	}
	return false;
}

bool Meshline::touches(Basisfunction *basis) const {
	if(span_u_line_) {
		if( (*basis)[1][0] < const_par_ && const_par_ < (*basis)[1][basis->getOrder(1)] &&
		   (start_ < (*basis)[0][basis->getOrder(0)]  || (*basis)[0][0] < stop_))
			return true;
	} else { // span-v line
		if( (*basis)[0][0] < const_par_ && const_par_ < (*basis)[0][basis->getOrder(0)] &&
		   (start_ < (*basis)[1][basis->getOrder(1)]  || (*basis)[1][0] < stop_))
			return true;
	}
	return false;
}

bool Meshline::splits(Basisfunction *basis) const {
	if(span_u_line_) {
		if( (*basis)[1][0] < const_par_ && const_par_ < (*basis)[1][basis->getOrder(1)] &&
		    start_-DOUBLE_TOL <= (*basis)[0][0]  && (*basis)[0][basis->getOrder(0)] <= stop_+DOUBLE_TOL)
			return  true;
	} else { // span-v line
		if( (*basis)[0][0] < const_par_ && const_par_ < (*basis)[0][basis->getOrder(0)] &&
		    start_-DOUBLE_TOL <= (*basis)[1][0]  && (*basis)[1][basis->getOrder(1)] <= stop_+DOUBLE_TOL)
			return  true;
	}
	return false;
}

/************************************************************************************************************************//**
 * \brief Checks if another meshline intersects this one
 * \param other  Pointer to the other meshline
 * \param parval Pointer (or null) to the place to store the parametric intersection value on *this line
 * \returns if *other intersects this line; touches at the start- or endpoints will return true
 ***************************************************************************************************************************/
bool Meshline::intersects(Meshline *other, double *parval) const {
	if(other->span_u_line_ == this->span_u_line_)
		return false;
	if(other->start_ <= this->const_par_  && this->const_par_  <= other->stop_ && 
	   this->start_  <= other->const_par_ && other->const_par_ <= this->stop_) {
		if(parval != nullptr)
			*parval = other->const_par_;
		return true;
	}
	return false;
}

bool Meshline::is_spanning_u() const {
	return span_u_line_;
}

bool Meshline::operator==(const Meshline &other) const {
	return span_u_line_  == other.span_u_line_ &&
	       const_par_    == other.const_par_ &&
	       start_        == other.start_ &&
	       stop_         == other.stop_ &&
	       multiplicity_ == other.multiplicity_;

}

// convenience macro for reading formated input
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing meshline\n"; exit(325); } ws(is); }
void Meshline::read(std::istream &is) {
	char nextChar;
	ws(is);
	nextChar = is.peek();
	if(nextChar == '[') { // first parametric direction interval => const v
		ASSERT_NEXT_CHAR('[');
		span_u_line_ = true;
		is >> start_;
		ASSERT_NEXT_CHAR(',');
		is >> stop_;
		ASSERT_NEXT_CHAR(']');
		ASSERT_NEXT_CHAR('x');
		is >> const_par_;
		ASSERT_NEXT_CHAR('(');
		is >> multiplicity_;
		ASSERT_NEXT_CHAR(')');
	} else {
		span_u_line_ = false;
		is >> const_par_;
		ASSERT_NEXT_CHAR('x');
		ASSERT_NEXT_CHAR('[');
		is >> start_;
		ASSERT_NEXT_CHAR(',');
		is >> stop_;
		ASSERT_NEXT_CHAR(']');
		ASSERT_NEXT_CHAR('(');
		is >> multiplicity_;
		ASSERT_NEXT_CHAR(')');
	}
}
#undef ASSERT_NEXT_CHAR

void Meshline::write(std::ostream &os) const {
	if(span_u_line_)
		os <<  "[" << start_ << ", " << stop_ << "] x " << const_par_ << " (" << multiplicity_ << ")";
	else // span-v line
		os << const_par_ << " x [" << start_ << ", " << stop_ << "] (" << multiplicity_ << ")";
}

void Meshline::writeMore(std::ostream &os) const {
	if(span_u_line_)
		os <<  "[" << start_ << ", " << stop_ << "] x " << const_par_ << " (" << multiplicity_ << ")";
	else // span-v line
		os << const_par_ << " x [" << start_ << ", " << stop_ << "] (" << multiplicity_ << ")";
	if(type_ == INITIAL)
		os << " INITIAL";
	else if(type_ == NEWLINE)
		os << " NEWLINE";
	else if(type_ == MERGING)
		os << " MERGING";
	else if(type_ == ELONGATION)
		os << " ELONGATION";
}

#undef DOUBLE_TOL

} // end namespace LR


