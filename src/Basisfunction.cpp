#include "LRSpline/Basisfunction.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Meshline.h"
#include <algorithm>
#include <cfloat>
#include <climits>
#include <set>

namespace LR {

Basisfunction::Basisfunction(int dim, int order_u, int order_v) {
	weight_       = 1;
	id_           = -1;
	knots_.resize(2);
	knots_[0].resize(order_u+1);
	knots_[1].resize(order_v+1);
	controlpoint_.resize(dim);
}

// Basisfunction::Basisfunction(RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 controlpoint, int dim, int order_u, int order_v, double weight) {
// }

Basisfunction::~Basisfunction() {
	PROFILE("Function destruction");
	for(uint i=0; i<support_.size(); i++)
		support_[i]->removeSupportFunction( (Basisfunction*) this);
}

Go::Point Basisfunction::getGrevilleParameter() const {
	Go::Point ans(2);
	ans *= 0;
	for(uint i=1; i<knots_[0].size()-1; i++)
		ans[0] += knots_[0][i];
	for(uint i=1; i<knots_[1].size()-1; i++)
		ans[1] += knots_[1][i];
	ans[0] /= (knots_[0].size()-2);
	ans[1] /= (knots_[1].size()-2);
	return ans;
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
	if(knots_[0][0] > u || u > knots_[0].back())
		return ;
	if(knots_[1][0] > v || v > knots_[1].back())
		return ;

	double ans_u[knots_[0].size()-1];
	double ans_v[knots_[1].size()-1];
	// double diff_u = 0;
	// double diff_v = 0;
	// double diff_2u[3];
	// double diff_2v[3];
	double **diff_u;
	double **diff_v;
	diff_u = new double*[derivs+1];
	diff_v = new double*[derivs+1];
	int diff_level;

	for(uint i=0; i<knots_[0].size()-1; i++) {
		if(u_from_right)
			ans_u[i] = (knots_[0][i] <= u && u <  knots_[0][i+1]) ? 1 : 0;
		else 
			ans_u[i] = (knots_[0][i] <  u && u <= knots_[0][i+1]) ? 1 : 0;
	}

	diff_level = knots_[0].size()-2;
	for(uint n=1; n<knots_[0].size()-1; n++, diff_level--) {
		if(diff_level <= derivs) {
			diff_u[diff_level] = new double[diff_level+1];
			for(int j=0; j<=diff_level; j++)
				diff_u[diff_level][j] = ans_u[j];
		}
		for(int d = diff_level; d<= derivs; d++) {
			for(uint j=0; j<knots_[0].size()-1-n; j++) {
				diff_u[d][j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (   n   )/(knots_[0][j+n]  -knots_[0][ j ])*diff_u[d][ j ];
				diff_u[d][j] -= (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (   n   )/(knots_[0][j+n+1]-knots_[0][j+1])*diff_u[d][j+1];
			}
		}

#if 0
// kept this to maybe make it easier to see the logic behind evaluating arbitrary high derivatives
		if(n==knots_[0].size()-2)
			for(int j=0; j<=knots_[0].size()-n; j++)
				diff_2u[j] = ans_u[j];
		if(n>=knots_[0].size()-2) {
			for(int j=0; j<knots_[0].size()-n; j++) {
				diff_2u[j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (   n   )/(knots_[0][j+n]  -knots_[0][ j ])*diff_2u[ j ];
				diff_2u[j] -= (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (   n   )/(knots_[0][j+n+1]-knots_[0][j+1])*diff_2u[j+1];
			}
		}
		if(n==knots_[0].size()-1) {
			int j=0;
			diff_u  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (   knots_[0].size()-1   )/(knots_[0][j+n]  -knots_[0][ j ])*ans_u[ j ];
			diff_u -= (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (   knots_[0].size()-1   )/(knots_[0][j+n+1]-knots_[0][j+1])*ans_u[j+1];
		}
#endif
		for(uint j=0; j<knots_[0].size()-1-n; j++) {
			ans_u[j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (  u-knots_[0][j]  )/(knots_[0][j+n]  -knots_[0][ j ])*ans_u[ j ];
			ans_u[j] += (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (knots_[0][j+n+1]-u)/(knots_[0][j+n+1]-knots_[0][j+1])*ans_u[j+1];
		}
	}
					   
	for(uint i=0; i<knots_[1].size()-1; i++) {
		if(v_from_right)
			ans_v[i] = (knots_[1][i] <= v && v <  knots_[1][i+1]) ? 1 : 0;
		else 
			ans_v[i] = (knots_[1][i] <  v && v <= knots_[1][i+1]) ? 1 : 0;
	}

	diff_level = knots_[1].size()-2;
	for(uint n=1; n<knots_[1].size()-1; n++, diff_level--) {
		if(diff_level <= derivs) {
			diff_v[diff_level] = new double[diff_level+1];
			for(int j=0; j<=diff_level; j++)
				diff_v[diff_level][j] = ans_v[j];
		}
		for(int d = diff_level; d<= derivs; d++) {
			for(uint j=0; j<knots_[1].size()-1-n; j++) {
				diff_v[d][j]  = (knots_[1][ j+n ]==knots_[1][ j ]) ? 0 : (   n   )/(knots_[1][j+n]  -knots_[1][ j ])*diff_v[d][ j ];
				diff_v[d][j] -= (knots_[1][j+n+1]==knots_[1][j+1]) ? 0 : (   n   )/(knots_[1][j+n+1]-knots_[1][j+1])*diff_v[d][j+1];
			}
		}
	
		for(uint j=0; j<knots_[1].size()-1-n; j++) {
			ans_v[j]  = (knots_[1][ j+n ]==knots_[1][ j ]) ? 0 : (  v-knots_[1][j]  )/(knots_[1][j+n]  -knots_[1][ j ])*ans_v[ j ];
			ans_v[j] += (knots_[1][j+n+1]==knots_[1][j+1]) ? 0 : (knots_[1][j+n+1]-v)/(knots_[1][j+n+1]-knots_[1][j+1])*ans_v[j+1];
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
	if(knots_[0][0] > u || u > knots_[0].back())
		return 0;
	if(knots_[1][0] > v || v > knots_[1].back())
		return 0;

	double ans_u[knots_[0].size()-1];
	double ans_v[knots_[1].size()-1];

	for(uint i=0; i<knots_[0].size()-1; i++) {
		if(u_from_right)
			ans_u[i] = (knots_[0][i] <= u && u <  knots_[0][i+1]) ? 1 : 0;
		else 
			ans_u[i] = (knots_[0][i] <  u && u <= knots_[0][i+1]) ? 1 : 0;
	}
	for(uint n=1; n<knots_[0].size()-1; n++)
		for(uint j=0; j<knots_[0].size()-1-n; j++) {
			ans_u[j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (  u-knots_[0][j]  )/(knots_[0][j+n]  -knots_[0][ j ])*ans_u[ j ];
			ans_u[j] += (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (knots_[0][j+n+1]-u)/(knots_[0][j+n+1]-knots_[0][j+1])*ans_u[j+1];
	}
					   
	for(uint i=0; i<knots_[1].size()-1; i++) {
		if(v_from_right)
			ans_v[i] = (knots_[1][i] <= v && v <  knots_[1][i+1]) ? 1 : 0;
		else 
			ans_v[i] = (knots_[1][i] <  v && v <= knots_[1][i+1]) ? 1 : 0;
	}
	for(uint n=1; n<knots_[1].size()-1; n++)
		for(uint j=0; j<knots_[1].size()-1-n; j++) {
			ans_v[j]  = (knots_[1][ j+n ]==knots_[1][ j ]) ? 0 : (  v-knots_[1][j]  )/(knots_[1][j+n]  -knots_[1][ j ])*ans_v[ j ];
			ans_v[j] += (knots_[1][j+n+1]==knots_[1][j+1]) ? 0 : (knots_[1][j+n+1]-v)/(knots_[1][j+n+1]-knots_[1][j+1])*ans_v[j+1];
	}

	return ans_u[0]*ans_v[0]*weight_;
}

void Basisfunction::getControlPoint(Go::Point &pt) const {
	pt.resize(controlpoint_.size());
	for(uint d=0; d<controlpoint_.size(); d++)
		pt[d] = controlpoint_[d];
}

bool Basisfunction::removeSupport(Element *el) {
	for(uint i=0; i<support_.size(); i++) {
		if(el == support_[i]) {
			support_[i] = support_.back();
			support_.pop_back();
			break;
		}
	}
	return false;
}

bool Basisfunction::addSupport(Element *el) {
	if(overlaps(el)) {
		support_.push_back(el);
		return true;
	}
	return false;
}

bool Basisfunction::overlaps(Element *el) const {
	return knots_[0][0]     < el->umax() &&
	       knots_[0].back() > el->umin() && 
	       knots_[1][0]     < el->vmax() &&
	       knots_[1].back() > el->vmin();
}

std::vector<Element*> Basisfunction::getExtendedSupport() {
	std::set<Element*> ans ;

	for(Element* e : support_)
		for(Basisfunction* b : e->support())
			for(Element* e2 : b->support())
				ans.insert(e2);
	std::vector<Element*> ans_vector(ans.size());
	std::copy(ans.begin(), ans.end(), ans_vector.begin());
	return ans_vector;
}

std::vector<Element*> Basisfunction::getMinimalExtendedSupport() {
	double min_du = DBL_MAX;
	double min_dv = DBL_MAX;
	Basisfunction *smallestGuy = NULL;

	bool edgeUmin = (knots_[0][0] == knots_[0][knots_[0].size()-2]);
	bool edgeUmax = (knots_[0][1] == knots_[0][knots_[0].size()-1]);
	bool edgeVmin = (knots_[1][0] == knots_[1][knots_[1].size()-2]);
	bool edgeVmax = (knots_[1][1] == knots_[1][knots_[1].size()-1]);

	if(! (edgeUmin || edgeUmax) )
		min_du = umax() - umin();
	if(! (edgeVmin || edgeVmax) )
		min_dv = vmax() - vmin();

	for(Element* e : support_) {
		for(Basisfunction* b : e->support()) {
			if( ((b->umin() <  this->umin() && b->umax() >= this->umax()) || 
			     (b->umin() <= this->umin() && b->umax() >  this->umax())) &&  // extend left OR right of the u-directions
			    ((b->vmin() <  this->vmin() && b->vmax() >= this->vmax()) || 
			     (b->vmin() <= this->vmin() && b->vmax() >  this->vmax())) &&  // AND extend up OR down in the v-direction
			    min_du >= b->umax()-b->umin() &&
			    min_dv >= b->vmax()-b->vmin()  ) { // AND less support than the last best guy

				min_du = b->umax()-b->umin();
				min_dv = b->vmax()-b->vmin();
				smallestGuy = b;
			}
		}
	}
	if(smallestGuy == NULL) {
		for(Element* e : support_) {
			for(Basisfunction* b : e->support()) {
				if(  b->umin() <= this->umin()    && 
				     b->umax() >= this->umax()    && 
				     b->vmin() <= this->vmin()    && 
				     b->vmax() >= this->vmax()    &&  // Cover *this' support
				    (b->umin() <  this->umin() || 
				     b->umax() >  this->umax() || 
				     b->vmin() <  this->vmin() ||    // Extend at least one direction
				     b->vmax() >  this->vmax() )  && 
					min_du >= b->umax()-b->umin() &&
					min_dv >= b->vmax()-b->vmin()  ) { // smaller support than the last best option

					min_du = b->umax()-b->umin();
					min_dv = b->vmax()-b->vmin();
					smallestGuy = b;
				}
			}
		}
	}

	std::vector<Element*> results(smallestGuy->nSupportedElements());
	std::copy(smallestGuy->supportedElementBegin(), smallestGuy->supportedElementEnd(), results.begin());
	return results;
}

void Basisfunction::setDimension(int dim) {
	controlpoint_.resize(dim);
	for(int i=0; i<dim; i++)
		controlpoint_[i]  = 0.0;
}

bool Basisfunction::operator==(const Basisfunction &other) const {
	for(uint i=0; i<=knots_[0].size()-1; i++)
		if(knots_[0][i] != other.knots_[0][i])
			return false;
	for(uint i=0; i<=knots_[1].size()-1; i++)
		if(knots_[1][i] != other.knots_[1][i])
			return false;
	return true;
}

void Basisfunction::operator+=(const Basisfunction &other) {
	double newWeight = weight_ + other.weight_;
	for(uint i=0; i<controlpoint_.size(); i++)
		controlpoint_[i] = (controlpoint_[i]*weight_ + other.controlpoint_[i]*other.weight_)/newWeight;
	weight_ = newWeight;
}

Basisfunction* Basisfunction::copy() {
    // Basisfunction *returnvalue = new Basisfunction(this->knots_[0], this->knots_[1], this->controlpoint_, controlpoint_.size(), knots_[0].size()+1, knots_[1].size()+1, this->weight_);
   
    // returnvalue->id_    = this->id_;

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
    returnvalue->knots_[0].size() = this->knots_[0].size();
    returnvalue->knots_[1].size() = this->knots_[1].size();
    returnvalue->weight_ = this->weight_;
    returnvalue->knots_[0] = this->knots_[0];
    returnvalue->knots_[1] = this->knots_[1];
    returnvalue->controlpoint_ = this->controlpoint_;
    */
    //std::vector<double> knots_;
    return NULL;
    
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
	for(uint i=0; i<knots_[0].size(); i++)
		is >> knots_[0][i];
	ASSERT_NEXT_CHAR(']');
	ASSERT_NEXT_CHAR('x');
	ASSERT_NEXT_CHAR('[');
	for(uint i=0; i<knots_[1].size(); i++)
		is >> knots_[1][i];
	ASSERT_NEXT_CHAR(']');

	// read control point
	for(uint i=0; i<controlpoint_.size(); i++)
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
	for(uint i=0; i<knots_[0].size(); i++)
		os << knots_[0][i] << " ";
	os << "] x [";
	for(uint i=0; i<knots_[1].size(); i++)
		os << knots_[1][i] << " ";
	os << "] ";
	for(uint i=0; i<controlpoint_.size(); i++)
		os << controlpoint_[i] << " ";
	os << "(" << weight_ << ")";
}

bool Basisfunction::isOverloaded()  const {
	for(uint i=0; i<support_.size(); i++) 
		if(! support_[i]->isOverloaded() )
			return false;
	return true;
}

int Basisfunction::getOverloadCount() const {
	int ans = INT_MAX;
	for(uint i=0; i<support_.size(); i++) 
		ans = (ans < support_[i]->getOverloadCount()) ? ans : support_[i]->getOverloadCount();
	return ans;
}

} // end namespace LR
