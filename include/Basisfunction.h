#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

namespace LR {

class Basisfunction {
public:
	Basisfunction(double *knot_u, double *knot_v, double *controlpoint, int dim, int order_u, int order_v, double weight=1.0);
	~Basisfunction();
	double evaluate(double u, double v) const;

private:
	double *knot_u_;
	double *knot_v_;
	double *controlpoint_;
	int dim_;
	int order_u_;
	int order_v_;
	double weight_;

};

} // end namespace LR

#endif

