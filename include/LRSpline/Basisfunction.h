#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

#include <GoTools/utils/Point.h>
#include <GoTools/geometry/Streamable.h>
#include <vector>

namespace LR {

// forms a 4bit binary mask where each bit corresponds to a side and corners have 2bits 'on'
enum parameterEdge {
NONE       = 0,
WEST       = 1,    // 0001
EAST       = 2,    // 0010
SOUTH      = 4,    // 0100
NORTH      = 8,    // 1000
SOUTH_WEST = 5,    // 0101
SOUTH_EAST = 7,    // 0110
NORTH_WEST = 9,    // 1001
NORTH_EAST = 10};  // 1010

typedef enum parameterEdge parameterEdge;

class Element;
class Meshline;

class Basisfunction : public Go::Streamable {
public:
	// constructors
	Basisfunction(int dim, int order_u, int order_v);
	template <typename RandomIterator1,
	          typename RandomIterator2,
	          typename RandomIterator3>
	Basisfunction(RandomIterator1 knot_u, RandomIterator2 knot_v, RandomIterator3 controlpoint, int dim, int order_u, int order_v, double weight=1.0) {
		weight_       = weight ;
		id_           = -1;
		knots_.resize(2);
		knots_[0].resize(order_u+1);
		knots_[1].resize(order_v+1);
		controlpoint_.resize(dim);

		std::copy(knot_u,       knot_u       + order_u+1,   knots_[0].begin());
		std::copy(knot_v,       knot_v       + order_v+1,   knots_[1].begin());
		std::copy(controlpoint, controlpoint + dim,         controlpoint_.begin());
	}
	~Basisfunction();
	Basisfunction* copy();

	double evaluate(double u, double v, bool u_from_right=true, bool v_from_right=true) const;
	void evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right=true, bool v_from_right=true) const;

	bool operator==(const Basisfunction &other) const;
	void operator+=(const Basisfunction &other) ;
	std::vector<double>& operator[](int i) { return knots_[i]; } ;

	bool overlaps(Element *el) const;
	bool addSupport(Element *el) ;
	bool removeSupport(Element *el) ;
	std::vector<Element*>::iterator supportedElementBegin(){ return support_.begin(); };
	std::vector<Element*>::iterator supportedElementEnd()  { return support_.end();   };
	const std::vector<Element*>     support() const        { return support_;         };

	std::vector<Element*> getExtendedSupport() ;
	std::vector<Element*> getMinimalExtendedSupport();

	bool isOverloaded() const;
	int getOverloadCount() const;

	// get/set methods
	void addPartialLine(Meshline *line);
	void getControlPoint(Go::Point &pt) const;
	void setId(int id)  { this->id_ = id; };
	void setDimension(int dim)  ;
	int getId() const   { return id_; };
	int nSupportedElements() { return support_.size(); };
	double dim()  const { return controlpoint_.size(); };
	double umin() const { return knots_[0][0];         };
	double umax() const { return knots_[0].back();     };
	double vmin() const { return knots_[1][0];         };
	double vmax() const { return knots_[1].back();     };
	int order_u() const { return knots_[0].size()-1;   };
	int order_v() const { return knots_[1].size()-1;   };
	std::vector<double>& getknots_u()              {return knots_[0]; };
	std::vector<double>& getknots_v()              {return knots_[1]; };
	std::vector<double>::iterator cp()             {return controlpoint_.begin(); };
	std::vector<double>::const_iterator cp() const {return controlpoint_.begin(); };
	double cp(int i)                         const {return controlpoint_[i]; };
	double w()                               const {return weight_; };
	Go::Point getGrevilleParameter() const;

	double grevilleParameter(int index, int order, std::vector<double> knot) const;
	void test(int index, int order, std::vector<double> knot){std::cout << "test" <<std::endl;};

	// IO-functions
	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

private:

	int    id_;
	double weight_;
	std::vector<double>               controlpoint_;
	std::vector<std::vector<double> > knots_;
	std::vector<Element*>             support_;

};

} // end namespace LR

#endif

