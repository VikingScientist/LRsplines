
#include "LRSplineSurface.h"

namespace LR {


LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational) {

	rational_ = rational;
	dim_      = dim;
	order_u_  = order_u;
	order_v_  = order_v;

	// basis_;
	// meshline_;
	// element_;
}

} // end namespace LR

