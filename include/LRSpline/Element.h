#ifndef ELEMENT_H
#define ELEMENT_H

#include <GoTools/geometry/Streamable.h>
#include <vector>
#include "HashSet.h"

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
	HashSet_iterator<Basisfunction*> supportBegin()             { return support_.begin(); };
	HashSet_iterator<Basisfunction*> supportEnd()               { return support_.end();   };
	HashSet_const_iterator<Basisfunction*> supportBegin() const { return support_.begin(); };
	HashSet_const_iterator<Basisfunction*> supportEnd()   const { return support_.end();   };
	const HashSet<Basisfunction*> support()               const { return support_;   };
	// Basisfunction* supportFunction(int i) { return support_[i];   };
	int nBasisFunctions() const           { return support_.size(); };
	void setId(int id)                    { this->id_ = id; };
	int getId() const                     { return id_; };
	void setUmin(double u)                           { start_u_ = u; };
	void setVmin(double v)                           { start_v_ = v; };
	void setUmax(double u)                           { stop_u_  = u; };
	void setVmax(double v)                           { stop_v_  = v; };

	bool isOverloaded() const;
	// int overloadedBasisCount() const;
	void resetOverloadCount()    { overloadCount = 0;      }
	int incrementOverloadCount() { return overloadCount++; }
	int getOverloadCount() const { return overloadCount;   }

	void updateBasisPointers(HashSet<Basisfunction*> &basis) ;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

private:
	double start_u_;
	double start_v_;
	double stop_u_;
	double stop_v_;
	int id_;

	HashSet<Basisfunction*> support_;
	std::vector<int> support_ids_; // temporary storage for the read() method only

	int overloadCount ;
	
};

} // end namespace LR

#endif

