#include "Basisfunction.h"
#include <algorithm>

namespace LR {

Basisfunction::Basisfunction(double *knot_u, double *knot_v, double *controlpoint, int dim, int order_u, int order_v, double weight) {
	dim_          = dim    ;
	order_u_      = order_u;
	order_v_      = order_v;
	weight_       = weight ;
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

double Basisfunction::evaluate(double u, double v) const {
	if(knot_u_[0] > u || u > knot_u_[order_u_])
		return 0;
	if(knot_v_[0] > v || v > knot_v_[order_v_])
		return 0;

	double ans_u[order_u_+1];
	double ans_v[order_v_+1];

	for(int i=0; i<order_u_; i++)
		if(knot_u_[i] <= u && u < knot_u_[i+1])
			ans_u[i] = 1.0;
	for(int n=1; n<order_u_+1; n++)
		for(int j=0; j<order_u_+1-n; j++)
			ans_u[j] = (  u-knot_u_[j]  )/(knot_u_[j+n]  -knot_u_[j]  )*ans_u[j]  + 
			           (knot_u_[j+n+1]-u)/(knot_u_[j+n+1]-knot_u_[j+1])*ans_u[j+1];
					   
	for(int i=0; i<order_v_; i++)
		if(knot_v_[i] <= v && v < knot_v_[i+1])
			ans_v[i] = 1.0;
	for(int n=1; n<order_v_+1; n++)
		for(int j=0; j<order_v_+1-n; j++)
			ans_v[j] = (  v-knot_v_[j]  )/(knot_v_[j+n]-  knot_v_[j])*ans_v[j]    +
			           (knot_v_[j+n+1]-v)/(knot_v_[j+n+1]-knot_v_[j+1])*ans_v[j+1];

	return ans_u[0]*ans_v[0]*weight_;
}

} // end namespace LR
