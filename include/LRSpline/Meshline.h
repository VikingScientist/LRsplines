#ifndef MESHLINE_H
#define MESHLINE_H

#include "Streamable.h"
#include <vector>

namespace LR {

enum meshlineExtension {
	ELONGATION,
	MERGING,
	NEWLINE,
	INITIAL
};

class Basisfunction;
class Element;

class Meshline : public Streamable {

public:
	Meshline();
	Meshline(bool span_u_line, double const_par, double start, double stop, int multiplicity);
	virtual ~Meshline();
	Meshline* copy();

	int nKnotsIn(Basisfunction *basis) const;
	bool splits(Basisfunction *basis) const;
	bool touches(Basisfunction *basis) const;
	bool splits(Element *el) const;
	bool touches(Element *el) const;
	bool intersects(Meshline *other, double *parval=nullptr) const;

	bool is_spanning_u() const;
	int multiplicity() const { return multiplicity_; };

	bool operator==(const Meshline &other) const;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;
	virtual void writeMore(std::ostream &os) const;

// private:
	bool span_u_line_;
	double const_par_;
	double start_;
	double stop_;
	int multiplicity_;

// even more private (only used for linear independence testing)
	enum meshlineExtension type_;

};

} // end namespace LR

#endif

