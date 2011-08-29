#ifndef MESHLINE_H
#define MESHLINE_H

#include <GoTools/geometry/Streamable.h>

namespace LR {

class Basisfunction;

class Meshline : public Go::Streamable {

public:
	Meshline(bool span_u_line, double const_par, double start, double stop, int multiplicity);

	bool containedIn(Basisfunction *basis) const;
	bool splits(Basisfunction *basis) const;

	bool is_spanning_u() const;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

// private:
	bool span_u_line_;
	double const_par_;
	double start_;
	double stop_;
	int multiplicity_;

};

} // end namespace LR

#endif

