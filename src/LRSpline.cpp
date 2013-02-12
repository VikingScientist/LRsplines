#include "LRSpline/LRSpline.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/Element.h"

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

} // end namespace LR
