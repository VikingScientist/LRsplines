#ifndef LR_SPLINE_H
#define LR_SPLINE_H

#include <vector>
#include <GoTools/utils/Point.h>
#include <GoTools/geometry/Streamable.h>
#include <GoTools/geometry/SplineSurface.h>


namespace LR {

class Basisfunction;
class Meshline;
class Element;

class LRSplineSurface : public Go::Streamable {

public:
	LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational=false);

	virtual void point(Go::Point &pt, double u, double v) const;
	virtual void point(std::vector<Go::Point> &pts, double upar, double vpar, int derivs) const;
	void computeBasis (double param_u, double param_v, Go::BasisPtsSf & result ) const;
	void computeBasis (double param_u, double param_v, Go::BasisDerivsSf & result ) const;

	void insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity=1);
	void insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity=1);
	
	void getGlobalKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const;

	virtual double startparam_u() const { return start_u_; };
	virtual double startparam_v() const { return start_v_; };
	virtual double endparam_u()   const { return end_u_; };
	virtual double endparam_v()   const { return end_v_; };

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;
	void writePostscriptMesh(std::ostream &out);

private:
	int split(bool insert_in_u, int function_index, double new_knot);
	void insert_line(bool const_u, double const_par, double start, double stop, int multiplicity);
	
	bool rational_;
	std::vector<Basisfunction*> basis_;
	std::vector<Meshline*> meshline_;
	std::vector<Element*> element_;
	int dim_;
	int order_u_;
	int order_v_;
	double start_u_;
	double start_v_;
	double end_u_;
	double end_v_;

};

} // end namespace LR

#endif

