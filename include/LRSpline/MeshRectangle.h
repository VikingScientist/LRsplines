#ifndef MESHRECTANGLE_H
#define MESHRECTANGLE_H

#include <GoTools/geometry/Streamable.h>
#include <vector>

namespace LR {

class Basisfunction;
class Element;

class MeshRectangle : public Go::Streamable {

public:
	MeshRectangle();
	MeshRectangle(double start_u,
	              double start_v,
	              double start_w,
	              double stop_u,
	              double stop_v,
	              double stop_w,
	              int multiplicity=1);
	~MeshRectangle();
	MeshRectangle* copy() const;

	int nKnotsIn(Basisfunction *basis) const;
	bool equals(const MeshRectangle *rect) const;
	bool overlaps(MeshRectangle *rect) const;
	bool splits(Basisfunction *basis) const;
	bool splits(Element *el) const;

	int    constDirection() const;
	double constParameter() const;

	bool operator==(const MeshRectangle &other) const;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

// private:
	std::vector<double> start_;
	std::vector<double> stop_;
	int multiplicity_;
	int constDir_;

};

} // end namespace LR

#endif

