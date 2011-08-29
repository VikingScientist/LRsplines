#include "Basisfunction.h"
#include <algorithm>

namespace LR {

Basisfunction::Basisfunction(double *knot_u, double *knot_v, double *controlpoint, int dim, int order_u, int order_v, double weight) {
	dim_          = dim    ;
	order_u_      = order_u;
	order_v_      = order_v;
	weight_       = weight ;
	edge_index_   = NONE;
	knot_u_       = new double[order_u+1];
	knot_v_       = new double[order_v+1];
	controlpoint_ = new double[dim];

	std::copy(knot_u,       knot_u       + order_u+1,   knot_u_);
	std::copy(knot_v,       knot_v       + order_v+1,   knot_v_);
	std::copy(controlpoint, controlpoint + dim,         controlpoint_);
}

Basisfunction::~Basisfunction() {
	delete[] knot_u_;
	delete[] knot_v_;
	delete[] controlpoint_;
}

void Basisfunction::evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right, bool v_from_right) const {
	if(derivs > 1) {
		std::cerr << "Basisfunction::evaluate() not implemented for more derivatives than 1\n";
		exit(325);
	}
	results.resize((derivs+1)*(derivs+2)/2);
	fill(results.begin(), results.end(), 0);
	if(knot_u_[0] > u || u > knot_u_[order_u_])
		return ;
	if(knot_v_[0] > v || v > knot_v_[order_v_])
		return ;

	double ans_u[order_u_];
	double ans_v[order_v_];
	double diff_u;
	double diff_v;

	for(int i=0; i<order_u_; i++) {
		if(u_from_right)
			ans_u[i] = (knot_u_[i] <= u && u <  knot_u_[i+1]) ? 1 : 0;
		else 
			ans_u[i] = (knot_u_[i] <  u && u <= knot_u_[i+1]) ? 1 : 0;
	}

	for(int n=1; n<order_u_; n++) {
		if(n==order_u_-1) {
			int j=0;
			diff_u  = (knot_u_[ j+n ]==knot_u_[ j ]) ? 0 : (   order_u_-1   )/(knot_u_[j+n]  -knot_u_[ j ])*ans_u[ j ];
			diff_u -= (knot_u_[j+n+1]==knot_u_[j+1]) ? 0 : (   order_u_-1   )/(knot_u_[j+n+1]-knot_u_[j+1])*ans_u[j+1];
		}
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
	for(int n=1; n<order_v_; n++) {
		if(n==order_v_-1) {
			int j=0;
			diff_v  = (knot_v_[ j+n ]==knot_v_[ j ]) ? 0 : (    order_v_-1  )/(knot_v_[j+n]  -knot_v_[ j ])*ans_v[ j ];
			diff_v -= (knot_v_[j+n+1]==knot_v_[j+1]) ? 0 : (    order_v_-1  )/(knot_v_[j+n+1]-knot_v_[j+1])*ans_v[j+1];
		}
		for(int j=0; j<order_v_-n; j++) {
			ans_v[j]  = (knot_v_[ j+n ]==knot_v_[ j ]) ? 0 : (  v-knot_v_[j]  )/(knot_v_[j+n]  -knot_v_[ j ])*ans_v[ j ];
			ans_v[j] += (knot_v_[j+n+1]==knot_v_[j+1]) ? 0 : (knot_v_[j+n+1]-v)/(knot_v_[j+n+1]-knot_v_[j+1])*ans_v[j+1];
		}
	}

	results[0] = ans_u[0] *ans_v[0] *weight_;
	if(derivs>0) {
		results[1] = diff_u  *ans_v[0]*weight_;
		results[2] = ans_u[0]*diff_v  *weight_;
	}
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

void Basisfunction::setEdge(parameterEdge edge_index) {
	edge_index_ = edge_index;
}

void Basisfunction::addEdge(parameterEdge edge_index) {
	edge_index_ = (parameterEdge) (edge_index_ | edge_index);
}

parameterEdge Basisfunction::getEdgeIndex() const {
	return edge_index_;
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

void Basisfunction::read(std::istream &is) {
}

void Basisfunction::write(std::ostream &os) const {
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
