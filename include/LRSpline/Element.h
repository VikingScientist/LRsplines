#ifndef ELEMENT_H
#define ELEMENT_H

#include <GoTools/geometry/Streamable.h>
#include <vector>

namespace LR {

class Basisfunction;
class Meshline;

class Element : public Go::Streamable {

public:
	Element();
	Element(double start_u, double start_v, double stop_u, double stop_v);
	void removeSupportFunction(Basisfunction *f);
	void addSupportFunction(Basisfunction *f);
	Element *split(bool split_u, double par_value);
	Element* copy();
	// get/set methods
	double umin() const         { return start_u_; };
	double vmin() const         { return start_v_; };
	double umax() const         { return stop_u_;  };
	double vmax() const         { return stop_v_;  };
	double area() const         { return (stop_v_-start_v_)*(stop_u_-start_u_);  };
	std::vector<Basisfunction*>::iterator supportBegin() { return support_.begin(); };
	std::vector<Basisfunction*>::iterator supportEnd()   { return support_.end();   };
	std::vector<Basisfunction*>::const_iterator supportBegin()const { return support_.begin(); };
	std::vector<Basisfunction*>::const_iterator supportEnd() const  { return support_.end();   };
	Basisfunction* supportFunction(int i) { return support_[i];   };
	int nBasisFunctions() const           { return support_.size(); };
	void setId(int id)                    { this->id_ = id; };
	int getId() const                     { return id_; };
	void setUmin(double u)                           { start_u_ = u; };
	void setVmin(double v)                           { start_v_ = v; };
	void setUmax(double u)                           { stop_u_  = u; };
	void setVmax(double v)                           { stop_v_  = v; };

	void addPartialLine(Meshline *line);
	void updateBasisPointers(std::vector<Basisfunction*> &basis) ;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

private:
	double start_u_;
	double start_v_;
	double stop_u_;
	double stop_v_;
	int id_;

	std::vector<Basisfunction*> support_;
	std::vector<int> support_ids_; // temporary storage for the read() method only
	
};

} // end namespace LR

#endif


