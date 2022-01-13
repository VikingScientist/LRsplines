
#include "LRSpline/MeshRectangle.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"
#include <algorithm>
#include <cmath>

namespace LR {

#define DOUBLE_TOL 1e-14

bool MeshRectangle::addUniqueRect(std::vector<MeshRectangle*> &rects, MeshRectangle* newRect) {
	for(MeshRectangle *m : rects) {
		if(m->equals(newRect) ) {
			return false;
		}
	}
	rects.push_back(newRect);
	return true;
};

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


MeshRectangle* MeshRectangle::copy() const {
	 MeshRectangle *returnvalue     = new MeshRectangle();

	 returnvalue->start_        = this->start_; // apperently vector::operator=() takes a deep copy
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
		if( fabs((*basis)[d][i] - start_[d]) < DOUBLE_TOL )
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

//! \brief returns true if \a rect overlaps or touches the edge of this MeshRectangle
bool MeshRectangle::overlaps(const MeshRectangle *rect) const {
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

//! \brief returns true if \a rect is completely contained in the range of this MeshRectangle
bool MeshRectangle::contains(const MeshRectangle *rect) const {
	int c1 = constDir_; // constant index
	int v1 = (c1+1)%3;  // first variable index
	int v2 = (c1+2)%3;  // second variable index

	if(constDir_ == rect->constDir_) {
		if(start_[c1] == rect->start_[c1] &&
		   start_[v1] <= rect->start_[v1] && stop_[v1] >= rect->stop_[v1] &&
		   start_[v2] <= rect->start_[v2] && stop_[v2] >= rect->stop_[v2] )
			return true;
	}

	return false;
}

int MeshRectangle::makeOverlappingRects(std::vector<MeshRectangle*> &newGuys, int meshIndex, bool allowSplits) {
	int c1  = constDir_; // constant index
	int v1  = (c1+1)%3;  // first variable index
	int v2  = (c1+2)%3;  // second variable index
	int v[] = {v1, v2};  // index way of referencing these
	bool addThisToNewGuys = false;
	MeshRectangle *rect = newGuys.at(meshIndex);
	if( ! this->overlaps(rect) )
		return 0;
	if( this->contains(rect) ) {
		newGuys.erase(newGuys.begin() + meshIndex);
		// std::cout << "Deleted: " << *rect << std::endl;
		// std::cout << "  contained in : " << *this << std::endl;
		delete rect;
		return 1;
	}
	if( rect->contains(this) ) {
		// std::cout << "Deleted: " << *this << std::endl;
		// std::cout << "  contained in : " << *rect << std::endl;
		return 2;
	}


	// for both variable indices... i=0 corresponds to checking everything left-right (direction v1),
	// i=1 is for checking everything up-down (direction v2).
	for(int i=0; i<2; i++) {
		int j = 1-i; // 0 or 1
		/*
		 *    MERGE RECTANGLES TO ONE
		 *  +-------+------+
		 *  |       |      |
		 *  |       |      |
		 *  |       |      |
		 *  +-------+------+
		 */
		// if the two mesh rectangles perfectly line up, keep only one of them
		if((stop_[v[i]]  == rect->start_[v[i]] ||
		    start_[v[i]] == rect->stop_[v[i]]  )) {

			if(start_[v[j]] == rect->start_[v[j]]  &&
			   stop_[v[j]]  == rect->stop_[v[j]]  ) {
				double min = std::min(start_[v[i]], rect->start_[v[i]]);
				double max = std::max(stop_[v[i]],  rect->stop_[v[i]] );
				// std::cout << "Deleted: " << *rect << std::endl;
				// std::cout << "  merged with  : " << *this << std::endl;
				newGuys.erase(newGuys.begin() + meshIndex);
				delete rect;
				start_[v[i]] = min;
				stop_[v[i]]  = max;
				if(addUniqueRect(newGuys, this))
					return 4;
				else
					return 5;
			/*
			 *    ADD ELONGATED RECT
			 * y3 +-------+
			 * y2 |       +------+         +-------+------+
			 *    |       |      |   =>    |              |
			 * y1 +-------+      |         +-------+------+
			 *            |      |                new one
			 * y0         +------+
			 *   x0      x1     x3
			 */

			} else if((start_[v[j]] < rect->start_[v[j]]   &&
			           stop_[v[j]]  < rect->stop_[v[j]] ) ||
			          (start_[v[j]] > rect->start_[v[j]]   &&
			           stop_[v[j]]  > rect->stop_[v[j]]  )) {
				double x0 = std::min(rect->start_[v[i]], start_[v[i]]);
				double x3 = std::max(rect->stop_[v[i]],  stop_[v[i]] );
				double y1 = std::max(rect->start_[v[j]], start_[v[j]]);
				double y2 = std::min(rect->stop_[v[j]],  stop_[v[j]] );
				double start[3];
				double stop[3];
				start[c1]   = start_[c1];
				stop[c1]    = stop_[c1];
				start[v[i]] = x0;
				stop[v[i]]  = x3;
				start[v[j]] = y1;
				stop[v[j]]  = y2;
				if(allowSplits) {
					MeshRectangle *m1 = new MeshRectangle(start, stop);
					if(!addUniqueRect(newGuys, m1))
						delete m1;
				}
			}
		}
		/*
		 *     EXPAND 'RECT' (RIGHT ONE)
		 *  +-------+
		 *  |    +--+------+
		 *  |    |  |      |
		 *  |    +--+------+
		 *  +-------+
		 */
		if(start_[v[i]] <  rect->start_[v[i]] &&
		   start_[v[j]] <= rect->start_[v[j]] &&
		   stop_[v[j]]  >= rect->stop_[v[j]]) { // expand the support of rect DOWN in v[i]
			rect->start_[v[i]] = start_[v[i]];
		}
		if(stop_[v[i]]  >  rect->stop_[v[i]]  &&
		   start_[v[j]] <= rect->start_[v[j]] &&
		   stop_[v[j]]  >= rect->stop_[v[j]]) { // expand the support of rect UP in v[i]
			rect->stop_[v[i]] = stop_[v[i]];
		}
		/*
		 *     EXPAND 'THIS' (LEFT ONE)
		 *        +------+
		 *  +-----+--+   |
		 *  |     |  |   |
		 *  +-----+--+   |
		 *        +------+
		 */
		if(rect->start_[v[i]] <  start_[v[i]] &&
		   rect->start_[v[j]] <= start_[v[j]] &&
		   rect->stop_[v[j]]  >= stop_[v[j]]) {
			start_[v[i]] = rect->start_[v[i]];
			addThisToNewGuys = true;
		}
		if(rect->stop_[v[i]]  >  stop_[v[i]]  &&
		   rect->start_[v[j]] <= start_[v[j]] &&
		   rect->stop_[v[j]]  >= stop_[v[j]]) {
			stop_[v[i]] = rect->stop_[v[i]];
			addThisToNewGuys = true;
		}
	}


	/*
	 *     DO NOTHING
	 *          +--+               +--+          +----+
	 *          |  |               |  |          |    |
	 *     +----+--+----+     +----+--+     +----+----+
	 *     |    |  |    |     |    |  |     |    |    |
	 *     +----+--+----+     +----+--+     |    |    |
	 *          |  |               |  |     +----+----+
	 *          +--+               +--+
	 *        note that this is after fixes above
	 */
	if((rect->stop_[v1]  >=       stop_[v1]   &&
	    rect->start_[v1] <=       start_[v1]  &&
	          stop_[v2]  >= rect->stop_[v2]   &&
	          start_[v2] <= rect->start_[v2]) ||
	   (      stop_[v1]  >= rect->stop_[v1]   &&
	          start_[v1] <= rect->start_[v1]  &&
	    rect->stop_[v2]  >=       stop_[v2]   &&
	    rect->start_[v2] <=       start_[v2])) {
		// ...
		;
	/*
	 *         MAKE TWO NEW ONES
	 *  y3  +-------+                        +--+
	 *      |       |                        |  |
	 *  y2  |    +--+----+           y2 +----+--+----+
	 *      |    |  |    |      =>      |    |  |    |
	 *  y1  +----+--+    |           y1 +----+--+----+
	 *           |       |                   |  |
	 *  y0       +-------+                   +--+
	 *     x0   x1 x2   x3                  x1  x2
	 *                                    these are the two new ones
	 */
	} else {
		double x0 = std::min(rect->start_[v1], start_[v1]);
		double x1 = std::max(rect->start_[v1], start_[v1]);
		double x2 = std::min(rect->stop_[v1],  stop_[v1] );
		double x3 = std::max(rect->stop_[v1],  stop_[v1] );
		double y0 = std::min(rect->start_[v2], start_[v2]);
		double y1 = std::max(rect->start_[v2], start_[v2]);
		double y2 = std::min(rect->stop_[v2],  stop_[v2] );
		double y3 = std::max(rect->stop_[v2],  stop_[v2] );
		double start[3];
		double stop[3];
		start[c1] = start_[c1];
		stop[c1]  = stop_[c1];

		if(allowSplits) {
			start[v1] = x0;
			stop[v1]  = x3;
			start[v2] = y1;
			stop[v2]  = y2;
			MeshRectangle *m1 = new MeshRectangle(start, stop);

			start[v1] = x1;
			stop[v1]  = x2;
			start[v2] = y0;
			stop[v2]  = y3;
			MeshRectangle *m2 = new MeshRectangle(start, stop);

			// std::cout << "Added: " << *m1 << std::endl;
			// std::cout << "Added: " << *m2 << std::endl;
			// std::cout << "  overlaps from  : " << *rect << std::endl;
			// std::cout << "  overlaps from  : " << *this << std::endl;

			if(!addUniqueRect(newGuys, m1))
				delete m1;
			if(!addUniqueRect(newGuys, m2))
				delete m2;
		}

	}
	if(addThisToNewGuys) {
		// std::cout << "Moved: " << *this << std::endl;
		if(addUniqueRect(newGuys, this))
			return 3;
		else
			return 6;
	}
	return 0;
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
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing meshrectangle\n"; exit(325); } ws(is); }
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

#undef DOUBLE_TOL

} // end namespace LR


