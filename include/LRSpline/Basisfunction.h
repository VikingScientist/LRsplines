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
	Basisfunction(const double *knot_u, const double *knot_v, double *controlpoint, int dim, int order_u, int order_v, double weight=1.0);
	~Basisfunction();
	Basisfunction* copy();

	double evaluate(double u, double v, bool u_from_right=true, bool v_from_right=true) const;
	void evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right=true, bool v_from_right=true) const;

	bool operator==(const Basisfunction &other) const;
	void operator+=(const Basisfunction &other) ;

	bool overlaps(Element *el) const;
	bool addSupport(Element *el) ;
	bool removeSupport(Element *el) ;
	std::vector<Element*>::iterator supportedElementBegin() ;
	std::vector<Element*>::iterator supportedElementEnd() ;
	std::vector<Meshline*>::iterator partialLineBegin() ;
	std::vector<Meshline*>::iterator partialLineEnd() ;

	std::vector<Element*> getExtendedSupport() ;
	std::vector<Element*> getMinimalExtendedSupport();

	void inheritPartialLine(Basisfunction *f);
	void removePartialLine(Meshline *m);

	bool isOverloaded() const;
	int getOverloadCount() const;

	// get/set methods
	void addPartialLine(Meshline *line);
	void getControlPoint(Go::Point &pt) const;
	void setId(int id)  { this->id_ = id; };
	void setDimension(int dim)  ;
	int getId() const   { return id_; };
	int nSupportedElements() { return support_.size(); };
	int nPartialLines() { return partial_line_.size(); };
	double umin() const { return knot_u_[0];        };
	double umax() const { return knot_u_[order_u_]; };
	double vmin() const { return knot_v_[0];        };
	double vmax() const { return knot_v_[order_v_]; };
	int order_u() const { return order_u_; };
	double *getknots_u() const {return knot_u_; };
	int order_v() const { return order_v_; };
	double *getknots_v() const {return knot_v_; };
	Go::Point getGrevilleParameter() const;
        double grevilleParameter(int index, int order, std::vector<double> knot) const;
        void test(int index, int order, std::vector<double> knot){std::cout << "test" <<std::endl;};
	// IO-functions
	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

// private:

	double *knot_u_;
	double *knot_v_;
	double *controlpoint_;
	int dim_;
	int order_u_;
	int order_v_;
	double weight_;
	std::vector<Element*> support_;
	std::vector<Meshline*> partial_line_;
	int id_;
	int order_;
	std::vector<double> knots_;

};

} // end namespace LR

#endif

