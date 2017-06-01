#include "LRSpline/LRSpline.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/Element.h"

typedef unsigned int uint;

namespace LR {

LRSpline::LRSpline() {
	dim_      = 0;
	element_.resize(0);
}

void LRSpline::generateIDs() const {
	uint i=0;
	for(Basisfunction *b : basis_)
		b->setId(i++);
	for(i=0; i<element_.size(); i++)
		element_[i]->setId(i);
}

void LRSpline::getEdgeFunctions(std::vector<Basisfunction*> &edgeFunctions, parameterEdge edge, int depth) const {
	edgeFunctions.clear();
	bool trivariate = (**basis_.begin()).nVariate() == 3;
	for(Basisfunction *b : basis_) {
		bool ok = true;
		if( edge & WEST )
			if((*b)[0][order_[0]-depth] != start_[0])
				ok = false;
		if( edge & EAST )
			if((*b)[0][depth] != end_[0])
				ok = false;
		if( edge & SOUTH )
			if((*b)[1][order_[1]-depth] != start_[1])
				ok = false;
		if( edge & NORTH )
			if((*b)[1][depth] != end_[1])
				ok = false;
		if( trivariate && (edge & BOTTOM) )
			if((*b)[2][order_[2]-depth] != start_[2])
				ok = false;
		if( trivariate && (edge & TOP ) )
			if((*b)[2][depth] != end_[2])
				ok = false;

		if(ok)
			edgeFunctions.push_back(b);
	}
}

void LRSpline::getEdgeElements( std::vector<Element*> &edgeElements, parameterEdge edge ) const {
	edgeElements.clear();
	bool trivariate = (**basis_.begin()).nVariate() == 3;
	for(Element * el : element_) {
		bool ok = true;
		if( edge & WEST )
			if(el->getParmin(0) != start_[0])
				ok = false;
		if( edge & EAST )
			if(el->getParmax(0) != end_[0])
				ok = false;
		if( edge & SOUTH )
			if(el->getParmin(1) != start_[1])
				ok = false;
		if( edge & NORTH )
			if(el->getParmax(1) != end_[1])
				ok = false;
		if( trivariate && (edge & BOTTOM) )
			if(el->getParmin(2) != start_[2])
				ok = false;
		if( trivariate && (edge & TOP ) )
			if(el->getParmax(2) != end_[2])
				ok = false;

		if(ok)
			edgeElements.push_back(el);
	}
}

bool LRSpline::setControlPoints(const std::vector<double>& controlpoints) {
	if((int) controlpoints.size() != dim_*basis_.size())
		return false;

	std::vector<double>::const_iterator newCP = controlpoints.begin();

	HashSet_iterator<Basisfunction*> bit;
	for(bit=basis_.begin(); bit!=basis_.end(); ++bit) {
		std::vector<double>::iterator cp = (**bit).cp();
		for(int i=0; i<dim_; i++, cp++, newCP++)
			*cp = *newCP;
	}
	return true;
}

void LRSpline::rebuildDimension(int dimvalue) {
	for(Basisfunction *b : basis_)
		b->setDimension(dimvalue);
	dim_ = dimvalue;
}


} // end namespace LR
