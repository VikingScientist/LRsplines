#ifndef LR_SPLINE_H
#define LR_SPLINE_H

#include <vector>

namespace LR {

class Basisfunction;
class Meshline;
class Element;

class LRSplineSurface {

public:
	LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational=false);

private:
	bool rational_;
	std::vector<Basisfunction*> basis_;
	std::vector<Meshline*> meshline_;
	std::vector<Element*> element_;
	int dim_;
	int order_u_;
	int order_v_;

};

} // end namespace LR

#endif



