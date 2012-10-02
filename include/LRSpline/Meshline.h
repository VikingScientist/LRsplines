#ifndef MESHLINE_H
#define MESHLINE_H

#include <GoTools/geometry/Streamable.h>
#include <vector>

namespace LR {

class Basisfunction;
class Element;

class Meshline : public Go::Streamable {

public:
	Meshline();
	Meshline(bool span_u_line, double const_par, double start, double stop, int multiplicity);
	~Meshline();
	Meshline* copy();

	int nKnotsIn(Basisfunction *basis) const;
	bool splits(Basisfunction *basis) const;
	bool touches(Basisfunction *basis) const;
	bool splits(Element *el) const;
	bool touches(Element *el) const;

	void addPartialTouch(Basisfunction *basis);
	void removePartialTouch(Basisfunction *basis);

	bool is_spanning_u() const;

	bool operator==(const Meshline &other) const;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

// private:
	bool span_u_line_;
	double const_par_;
	double start_;
	double stop_;
	int multiplicity_;
	std::vector<Basisfunction*> touching_Basisfunction;

};

} // end namespace LR

#endif

