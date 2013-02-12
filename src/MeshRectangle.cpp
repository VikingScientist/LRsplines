
#include "LRSpline/MeshRectangle.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

namespace LR {

#define DOUBLE_TOL 1e-14
#define MY_STUPID_FABS(x) (((x)>0)?(x):-(x))

MeshRectangle::MeshRectangle() {
	start_.resize(3);
	stop_.resize(3);
	start_[0]     = 0;
	start_[1]     = 0;
	start_[2]     = 0;
	stop_[0]      = 0;
	stop_[1]      = 0;
	stop_[2]      = 0;
	multiplicity_ = 0;
	constDir_     = -1;
}

MeshRectangle::MeshRectangle(double start_u,
                             double start_v,
                             double start_w,
                             double stop_u,
                             double stop_v,
                             double stop_w,
                             int multiplicity) {
	start_.resize(3);
	stop_.resize(3);
	start_[0]     = start_u;
	start_[1]     = start_v;
	start_[2]     = start_w;
	stop_[0]      = stop_u; 
	stop_[1]      = stop_v; 
	stop_[2]      = stop_w; 
	multiplicity_ =  multiplicity;
	constDir_     = -1;

	for(int i=0; i<3; i++)
		if(start_[i] == stop_[i])
			constDir_ = i;
	if(constDir_ == -1)
		std::cerr << "Error creating MeshRectangle: Not parallel to the parametric axis\n";
}

MeshRectangle::~MeshRectangle() {
}


MeshRectangle* MeshRectangle::copy() {
	 MeshRectangle *returnvalue     = new MeshRectangle();
	 
	 returnvalue->start_        = this->start_;
	 returnvalue->stop_         = this->stop_ ;
	 returnvalue->multiplicity_ = this->multiplicity_;
	 returnvalue->constDir_     = this->constDir_;
	 return returnvalue;
}

int MeshRectangle::constDirection() const {
	return constDir_;
}

double MeshRectangle::constParameter() const {
	if(constDir_ > 2 || constDir_ < 0) return -1.0;
	return start_[constDir_];
}

int MeshRectangle::nKnotsIn(Basisfunction *basis) const {
	int hits = 0;
	int d    = constDir_;
	for(int i=0; i<=basis->getOrder(d); i++)
		if( MY_STUPID_FABS((*basis)[d][i] - start_[d]) < DOUBLE_TOL )
			hits++;
	return hits;
}

bool MeshRectangle::equals(const MeshRectangle *rect) const {
	for(int i=0; i<3; i++) {
		if(start_[i] != rect->start_[i])     return false;
		if(stop_[i]  != rect->stop_[i])      return false;
	}
	if(multiplicity_ != rect->multiplicity_) return false;

	return true;
}

bool MeshRectangle::overlaps(MeshRectangle *rect) const {
	int c1 = constDir_; // constant index
	int v1 = (c1+1)%3;  // first variable index
	int v2 = (c1+2)%3;  // second variable index

	if(constDir_ == rect->constDir_) {
		if(start_[c1] == rect->start_[c1] &&
		   start_[v1] <= rect->stop_[v1] && stop_[v1] >= rect->start_[v1] &&
		   start_[v2] <= rect->stop_[v2] && stop_[v2] >= rect->start_[v2])
			return true;
	}

	return false;
}

bool MeshRectangle::splits(Element *el) const {
	if(constDir_ == -1) {
		std::cerr << "MeshRectangle::splits(Element*) - constDir_ not set on meshrectangle\n";
		exit(3421);
	}
	
	int c1 = constDir_; // constant index
	int v1 = (c1+1)%3;  // first variable index
	int v2 = (c1+2)%3;  // second variable index

	if(el->getParmin(c1) < start_[c1]  && start_[c1] < el->getParmax(c1) && 
	   start_[v1] <= el->getParmin(v1) && el->getParmax(v1) <= stop_[v1] && 
	   start_[v2] <= el->getParmin(v2) && el->getParmax(v2) <= stop_[v2])
		return true;

	return false;
}

bool MeshRectangle::splits(Basisfunction *basis) const {
	if(constDir_ == -1) {
		std::cerr << "MeshRectangle::splits(Basisfunction*) - constDir_ not set on meshrectangle\n";
		exit(3421);
	}
	
	int c1 = constDir_; // constant index
	int v1 = (c1+1)%3;  // first variable index
	int v2 = (c1+2)%3;  // second variable index

	if(basis->getParmin(c1) < start_[c1]  && start_[c1] < basis->getParmax(c1) && 
	   start_[v1] <= basis->getParmin(v1) && basis->getParmax(v1) <= stop_[v1] && 
	   start_[v2] <= basis->getParmin(v2) && basis->getParmax(v2) <= stop_[v2])
		return true;

	return false;
}

bool MeshRectangle::operator==(const MeshRectangle &other) const {
	return this->equals(&other);
}

// convenience macro for reading formated input
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing meshrectangle\n"; std::cout << is; exit(325); } ws(is); }
void MeshRectangle::read(std::istream &is) {
	char nextChar;
	ws(is);

	ASSERT_NEXT_CHAR('[');
	is >> start_[0];
	ASSERT_NEXT_CHAR(',');
	is >> stop_[0];
	ASSERT_NEXT_CHAR(']');
	ASSERT_NEXT_CHAR('x');

	ASSERT_NEXT_CHAR('[');
	is >> start_[1];
	ASSERT_NEXT_CHAR(',');
	is >> stop_[1];
	ASSERT_NEXT_CHAR(']');
	ASSERT_NEXT_CHAR('x');

	ASSERT_NEXT_CHAR('[');
	is >> start_[2];
	ASSERT_NEXT_CHAR(',');
	is >> stop_[2];
	ASSERT_NEXT_CHAR(']');

	ASSERT_NEXT_CHAR('(');
	is >> multiplicity_;
	ASSERT_NEXT_CHAR(')');


	for(int i=0; i<3; i++)
		if(start_[i] == stop_[i])
			constDir_ = i;
	if(constDir_ == -1)
		std::cerr << "Error creating MeshRectangle: Not parallel to the parametric axis\n";
}
#undef ASSERT_NEXT_CHAR

void MeshRectangle::write(std::ostream &os) const {
	os << "[" << start_[0] << ", " << stop_[0] << "] x ";
	os << "[" << start_[1] << ", " << stop_[1] << "] x ";
	os << "[" << start_[2] << ", " << stop_[2] << "] ";
	os << "(" << multiplicity_ << ")";
}

#undef MY_STUPID_FABS
#undef DOUBLE_TOL

} // end namespace LR


