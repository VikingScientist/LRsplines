#include "LRSpline/Basisfunction.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Meshline.h"
#include <algorithm>

namespace LR {

Basisfunction::Basisfunction(int dim, int order_u, int order_v) {
	dim_          = dim    ;
	order_u_      = order_u;
	order_v_      = order_v;
	weight_       = 1;
	id_           = -1;
	knot_u_       = new double[order_u+1];
	knot_v_       = new double[order_v+1];
	controlpoint_ = new double[dim];
}

Basisfunction::Basisfunction(const double *knot_u, const double *knot_v, double *controlpoint, int dim, int order_u, int order_v, double weight) {
	dim_          = dim    ;
	order_u_      = order_u;
	order_v_      = order_v;
	weight_       = weight ;
	id_           = -1;
	knot_u_       = new double[order_u+1];
	knot_v_       = new double[order_v+1];
	controlpoint_ = new double[dim];

	std::copy(knot_u,       knot_u       + order_u+1,   knot_u_);
	std::copy(knot_v,       knot_v       + order_v+1,   knot_v_);
	std::copy(controlpoint, controlpoint + dim,         controlpoint_);
}

Basisfunction::~Basisfunction() {
	PROFILE("Function destruction");
	for(uint i=0; i<support_.size(); i++)
		support_[i]->removeSupportFunction( (Basisfunction*) this);
	// for(uint i=0; i<partial_line_.size(); i++)
		// partial_line_[i]->removePartialTouch( (Basisfunction*) this);
	delete[] knot_u_;
	delete[] knot_v_;
	delete[] controlpoint_;
}


double Basisfunction::grevilleParameter(int index, int order, std::vector<double> knot) const
{
   double greville = 0.0;
   for (int i = 1; i < order; ++i)
      greville += knot[index+i];

   return greville/(order - 1);
}


void Basisfunction::evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right, bool v_from_right) const {
	results.resize((derivs+1)*(derivs+2)/2);
	fill(results.begin(), results.end(), 0);
	if(knot_u_[0] > u || u > knot_u_[order_u_])
		return ;
	if(knot_v_[0] > v || v > knot_v_[order_v_])
		return ;

	double ans_u[order_u_];
	double ans_v[order_v_];
	// double diff_u = 0;
	// double diff_v = 0;
	// double diff_2u[3];
	// double diff_2v[3];
	double **diff_u;
	double **diff_v;
	diff_u = new double*[derivs+1];
	diff_v = new double*[derivs+1];
	int diff_level;

	for(int i=0; i<order_u_; i++) {
		if(u_from_right)
			ans_u[i] = (knot_u_[i] <= u && u <  knot_u_[i+1]) ? 1 : 0;
		else 
			ans_u[i] = (knot_u_[i] <  u && u <= knot_u_[i+1]) ? 1 : 0;
	}

	diff_level = order_u_-1;
	for(int n=1; n<order_u_; n++, diff_level--) {
		if(diff_level <= derivs) {
			diff_u[diff_level] = new double[diff_level+1];
			for(int j=0; j<=diff_level; j++)
				diff_u[diff_level][j] = ans_u[j];
		}
		for(int d = diff_level; d<= derivs; d++) {
			for(int j=0; j<order_u_-n; j++) {
				diff_u[d][j]  = (knot_u_[ j+n ]==knot_u_[ j ]) ? 0 : (   n   )/(knot_u_[j+n]  -knot_u_[ j ])*diff_u[d][ j ];
				diff_u[d][j] -= (knot_u_[j+n+1]==knot_u_[j+1]) ? 0 : (   n   )/(knot_u_[j+n+1]-knot_u_[j+1])*diff_u[d][j+1];
			}
		}

#if 0
// kept this to maybe make it easier to see the logic behind evaluating arbitrary high derivatives
		if(n==order_u_-2)
			for(int j=0; j<=order_u_-n; j++)
				diff_2u[j] = ans_u[j];
		if(n>=order_u_-2) {
			for(int j=0; j<order_u_-n; j++) {
				diff_2u[j]  = (knot_u_[ j+n ]==knot_u_[ j ]) ? 0 : (   n   )/(knot_u_[j+n]  -knot_u_[ j ])*diff_2u[ j ];
				diff_2u[j] -= (knot_u_[j+n+1]==knot_u_[j+1]) ? 0 : (   n   )/(knot_u_[j+n+1]-knot_u_[j+1])*diff_2u[j+1];
			}
		}
		if(n==order_u_-1) {
			int j=0;
			diff_u  = (knot_u_[ j+n ]==knot_u_[ j ]) ? 0 : (   order_u_-1   )/(knot_u_[j+n]  -knot_u_[ j ])*ans_u[ j ];
			diff_u -= (knot_u_[j+n+1]==knot_u_[j+1]) ? 0 : (   order_u_-1   )/(knot_u_[j+n+1]-knot_u_[j+1])*ans_u[j+1];
		}
#endif
		for(int j=0; j<order_u_-n; j++) {
			ans_u[j]  = (knot_u_[ j+n ]==knot_u_[ j ]) ? 0 : (  u-knot_u_[j]  )/(knot_u_[j+n]  -knot_u_[ j ])*ans_u[ j ];
			ans_u[j] += (knot_u_[j+n+1]==knot_u_[j+1]) ? 0 : (knot_u_[j+n+1]-u)/(knot_u_[j+n+1]-knot_u_[j+1])*ans_u[j+1];
		}
	}
					   
	for(int i=0; i<order_v_; i++) {
		if(v_from_right)
			ans_v[i] = (knot_v_[i] <= v && v <  knot_v_[i+1]) ? 1 : 0;
		else 
			ans_v[i] = (knot_v_[i] <  v && v <= knot_v_[i+1]) ? 1 : 0;
	}

	diff_level = order_v_-1;
	for(int n=1; n<order_v_; n++, diff_level--) {
		if(diff_level <= derivs) {
			diff_v[diff_level] = new double[diff_level+1];
			for(int j=0; j<=diff_level; j++)
				diff_v[diff_level][j] = ans_v[j];
		}
		for(int d = diff_level; d<= derivs; d++) {
			for(int j=0; j<order_v_-n; j++) {
				diff_v[d][j]  = (knot_v_[ j+n ]==knot_v_[ j ]) ? 0 : (   n   )/(knot_v_[j+n]  -knot_v_[ j ])*diff_v[d][ j ];
				diff_v[d][j] -= (knot_v_[j+n+1]==knot_v_[j+1]) ? 0 : (   n   )/(knot_v_[j+n+1]-knot_v_[j+1])*diff_v[d][j+1];
			}
		}
	
		for(int j=0; j<order_v_-n; j++) {
			ans_v[j]  = (knot_v_[ j+n ]==knot_v_[ j ]) ? 0 : (  v-knot_v_[j]  )/(knot_v_[j+n]  -knot_v_[ j ])*ans_v[ j ];
			ans_v[j] += (knot_v_[j+n+1]==knot_v_[j+1]) ? 0 : (knot_v_[j+n+1]-v)/(knot_v_[j+n+1]-knot_v_[j+1])*ans_v[j+1];
		}
	}

	// collect results
	diff_u[0] = ans_u;
	diff_v[0] = ans_v;
	int ip=0;
	for(int totDeriv=0; totDeriv<=derivs; totDeriv++)
		for(int vDeriv=0; vDeriv<=totDeriv; vDeriv++)
			results[ip++] = diff_u[totDeriv-vDeriv][0] * diff_v[vDeriv][0] * weight_;
	
	// clean up memory allocations
	for(int i=1; i<derivs+1; i++) {
		delete[] diff_u[i];
		delete[] diff_v[i];
	}
	delete[] diff_u;
	delete[] diff_v;
	
}

double Basisfunction::evaluate(double u, double v, bool u_from_right, bool v_from_right) const {
	if(knot_u_[0] > u || u > knot_u_[order_u_])
		return 0;
	if(knot_v_[0] > v || v > knot_v_[order_v_])
		return 0;

	double ans_u[order_u_];
	double ans_v[order_v_];

	for(int i=0; i<order_u_; i++) {
		if(u_from_right)
			ans_u[i] = (knot_u_[i] <= u && u <  knot_u_[i+1]) ? 1 : 0;
		else 
			ans_u[i] = (knot_u_[i] <  u && u <= knot_u_[i+1]) ? 1 : 0;
	}
	for(int n=1; n<order_u_; n++)
		for(int j=0; j<order_u_-n; j++) {
			ans_u[j]  = (knot_u_[ j+n ]==knot_u_[ j ]) ? 0 : (  u-knot_u_[j]  )/(knot_u_[j+n]  -knot_u_[ j ])*ans_u[ j ];
			ans_u[j] += (knot_u_[j+n+1]==knot_u_[j+1]) ? 0 : (knot_u_[j+n+1]-u)/(knot_u_[j+n+1]-knot_u_[j+1])*ans_u[j+1];
	}
					   
	for(int i=0; i<order_v_; i++) {
		if(v_from_right)
			ans_v[i] = (knot_v_[i] <= v && v <  knot_v_[i+1]) ? 1 : 0;
		else 
			ans_v[i] = (knot_v_[i] <  v && v <= knot_v_[i+1]) ? 1 : 0;
	}
	for(int n=1; n<order_v_; n++)
		for(int j=0; j<order_v_-n; j++) {
			ans_v[j]  = (knot_v_[ j+n ]==knot_v_[ j ]) ? 0 : (  v-knot_v_[j]  )/(knot_v_[j+n]  -knot_v_[ j ])*ans_v[ j ];
			ans_v[j] += (knot_v_[j+n+1]==knot_v_[j+1]) ? 0 : (knot_v_[j+n+1]-v)/(knot_v_[j+n+1]-knot_v_[j+1])*ans_v[j+1];
	}

	return ans_u[0]*ans_v[0]*weight_;
}

void Basisfunction::getControlPoint(Go::Point &pt) const {
	pt.resize(dim_);
	for(int d=0; d<dim_; d++)
		pt[d] = controlpoint_[d];
}

std::vector<Meshline*>::iterator Basisfunction::partialLineBegin()  {
	return partial_line_.begin();
}

std::vector<Meshline*>::iterator Basisfunction::partialLineEnd()  {
	return partial_line_.end();
}

std::vector<Element*>::iterator Basisfunction::supportedElementBegin()  {
	return support_.begin();
}

std::vector<Element*>::iterator Basisfunction::supportedElementEnd()  {
	return support_.end();
}

bool Basisfunction::removeSupport(Element *el) {
	for(uint i=0; i<support_.size(); i++) {
		if(el == support_[i]) {
			support_[i] = support_.back();
			support_.pop_back();
			break;
		}
	}
	// if(overlaps(el)) {
		// support_.push_back(el);
		// return true;
	// }
	return false;
}

bool Basisfunction::addSupport(Element *el) {
	if(overlaps(el)) {
		// for(uint i=0; i<support_.size(); i++)
			// if(support_[i] == el) return true;
		support_.push_back(el);
		return true;
	}
	return false;
}

bool Basisfunction::overlaps(Element *el) const {
	return knot_u_[0]        < el->umax() &&
	       knot_u_[order_u_] > el->umin() && 
	       knot_v_[0]        < el->vmax() &&
	       knot_v_[order_v_] > el->vmin();
}

void Basisfunction::inheritPartialLine(Basisfunction *f) {
	std::vector<Meshline*>::iterator it;
	for(it=f->partialLineBegin(); it != f->partialLineEnd(); it++) {
		if((*it)->touches(this))  {
			partial_line_.push_back(*it);
			// (*it)->addPartialTouch((Basisfunction*) this);
		}
	}
}

void Basisfunction::addPartialLine(Meshline *line) {
	// partial_line_.push_back(line);
	// line->addPartialTouch((Basisfunction*) this); 
};

void Basisfunction::removePartialLine(Meshline *m) {
	for(uint i=0; i<partial_line_.size(); i++) {
		if(*partial_line_[i] == *m) {
			partial_line_[i]->removePartialTouch( (Basisfunction*) this);
			partial_line_.erase(partial_line_.begin() + i);
			return;
		}
	}
}

bool Basisfunction::operator==(const Basisfunction &other) const {
	for(int i=0; i<=order_u_; i++)
		if(knot_u_[i] != other.knot_u_[i])
			return false;
	for(int i=0; i<=order_v_; i++)
		if(knot_v_[i] != other.knot_v_[i])
			return false;
	return true;
}

void Basisfunction::operator+=(const Basisfunction &other) {
	double newWeight = weight_ + other.weight_;
	for(int i=0; i<dim_; i++)
		controlpoint_[i] = (controlpoint_[i]*weight_ + other.controlpoint_[i]*other.weight_)/newWeight;
	weight_ = newWeight;
}

Basisfunction* Basisfunction::copy()
	  {
	    Basisfunction *returnvalue = new Basisfunction(this->knot_u_, this->knot_v_, this->controlpoint_,this->dim_,this->order_u_,this->order_v_,this->weight_);

	   
	    returnvalue->id_ = this->id_;
	    returnvalue->order_ = this->order_;
	    returnvalue->knots_ = this->knots_;
	    

	    /*
	    for(int i=0;i<support_.size(); i++)
	      {
		returnvalue -> support_.push_back(this->support_[i]->copy());
	      }
	    */

	    //loop
	    //returnvalue->knots_ = this->knots_;

	    //loop
	    /*
	    returnvalue->dim_ = this->dim_;
	    returnvalue->order_u_ = this->order_u_;
	    returnvalue->order_v_ = this->order_v_;
	    returnvalue->weight_ = this->weight_;
	    returnvalue->knot_u_ = this->knot_u_;
	    returnvalue->knot_v_ = this->knot_v_;
	    returnvalue->controlpoint_ = this->controlpoint_;
	    */
	    //std::vector<double> knots_;
	    return returnvalue;
	    
	  }








// convenience macro for reading formated input
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing basis function\n"; exit(324); } ws(is); }
void Basisfunction::read(std::istream &is) {
	char nextChar;

	// read id tag
	is >> id_;
	ws(is);
	ASSERT_NEXT_CHAR(':');

	// read knot vectors
	ASSERT_NEXT_CHAR('[');
	for(int i=0; i<=order_u_; i++)
		is >> knot_u_[i];
	ASSERT_NEXT_CHAR(']');
	ASSERT_NEXT_CHAR('x');
	ASSERT_NEXT_CHAR('[');
	for(int i=0; i<=order_v_; i++)
		is >> knot_v_[i];
	ASSERT_NEXT_CHAR(']');

	// read control point
	for(int i=0; i<dim_; i++)
		is >> controlpoint_[i];

	// read weight
	ASSERT_NEXT_CHAR('(');
	is >> weight_;
	ASSERT_NEXT_CHAR(')');
}
#undef ASSERT_NEXT_CHAR

void Basisfunction::write(std::ostream &os) const {
	os << id_ << ":";
	os << "[";
	for(int i=0; i<order_u_+1; i++)
		os << knot_u_[i] << " ";
	os << "] x [";
	for(int i=0; i<order_v_+1; i++)
		os << knot_v_[i] << " ";
	os << "] ";
	for(int i=0; i<dim_; i++)
		os << controlpoint_[i] << " ";
	os << "(" << weight_ << ")";
}


} // end namespace LR
