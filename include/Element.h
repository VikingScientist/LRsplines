#ifndef ELEMENT_H
#define ELEMENT_H

#include <GoTools/geometry/Streamable.h>
#include <vector>

namespace LR {

class Basisfunction;

class Element : public Go::Streamable {

public:
	Element(double start_u, double start_v, double stop_u, double stop_v);
	void removeSupportFunction(Basisfunction *f);
	void addSupportFunction(Basisfunction *f);
	Element *split(bool split_u, double par_value);

	double umin() const { return start_u_; };
	double vmin() const { return start_v_; };
	double umax() const { return stop_u_;  };
	double vmax() const { return stop_v_;  };
	std::vector<Basisfunction*>::iterator supportBegin() { return support_.begin(); };
	std::vector<Basisfunction*>::iterator supportEnd()   { return support_.end();   };

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

private:
	double start_u_;
	double start_v_;
	double stop_u_;
	double stop_v_;

	std::vector<Basisfunction*> support_;
	
};

} // end namespace LR

#endif


