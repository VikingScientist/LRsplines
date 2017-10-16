#ifndef MESHRECTANGLE_H
#define MESHRECTANGLE_H

#include "Streamable.h"
#include <vector>

namespace LR {

class Basisfunction;
class Element;

class MeshRectangle : public Streamable {

public:
	MeshRectangle();
	MeshRectangle(double start_u,
	              double start_v,
	              double start_w,
	              double stop_u,
	              double stop_v,
	              double stop_w,
	              int multiplicity=1);

	template <typename RandomIterator1,
	          typename RandomIterator2>
	MeshRectangle(RandomIterator1 start, RandomIterator2 stop, int multiplicity=1) {
		multiplicity_ = multiplicity;
		start_.resize(3);
		stop_.resize(3);
		std::copy(start, start+3, start_.begin());
		std::copy(stop,  stop+3,  stop_.begin() );
		for(int i=0; i<3; i++)
			if(start_[i] == stop_[i])
				constDir_ = i;
		if(constDir_ == -1)
			std::cerr << "Error creating MeshRectangle: Not parallel to the parametric axis\n";
	}
	virtual ~MeshRectangle();
	MeshRectangle* copy() const;

	int nKnotsIn(Basisfunction *basis) const;
	bool equals(const MeshRectangle *rect) const;
	bool overlaps(const MeshRectangle *rect) const;
	bool contains(const MeshRectangle *rect) const;
	bool splits(Basisfunction *basis) const;
	bool splits(Element *el) const;
	int makeOverlappingRects(std::vector<MeshRectangle*> &newGuys, int meshIndex, bool allowSplits) ;

	int    multiplicity()   const { return multiplicity_; };
	int    constDirection() const;
	double constParameter() const;

	bool operator==(const MeshRectangle &other) const;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

	static bool addUniqueRect(std::vector<MeshRectangle*> &rects, MeshRectangle* newRect);

// private:
	std::vector<double> start_;
	std::vector<double> stop_;
	int multiplicity_;
	int constDir_;

};

} // end namespace LR

#endif

