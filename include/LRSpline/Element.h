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
	Element(int dim);
	Element(double start_u, double start_v, double stop_u, double stop_v);
	Element(std::vector<double> &lowerLeft, std::vector<double> &upperRight);
	void removeSupportFunction(Basisfunction *f);
	void addSupportFunction(Basisfunction *f);
	Element *split(int splitDim, double par_value);
	Element* copy();
	// get/set methods
	double umin() const         { return min[0]; };
	double vmin() const         { return min[1]; };
	double umax() const         { return max[0]; };
	double vmax() const         { return max[1]; };
	double area() const         { return (max[1]-min[1])*(max[0]-min[0]);  };
	HashSet_iterator<Basisfunction*> supportBegin()                 { return support_.begin(); };
	HashSet_iterator<Basisfunction*> supportEnd()                   { return support_.end();   };
	HashSet_const_iterator<Basisfunction*> constSupportBegin()const { return support_.begin(); };
	HashSet_const_iterator<Basisfunction*> constSupportEnd()  const { return support_.end();   };
	const HashSet<Basisfunction*>& support()                  const { return support_;   };
	// Basisfunction* supportFunction(int i) { return support_[i];   };
	int nBasisFunctions() const           { return support_.size(); };
	void setId(int id)                    { this->id_ = id; };
	int  getId() const                    { return id_; };
	int  getDim() const                   { return min.size(); };
	void setUmin(double u)                { min[0] = u; };
	void setVmin(double v)                { min[1] = v; };
	void setUmax(double u)                { max[0] = u; };
	void setVmax(double v)                { max[1] = v; };

	bool isOverloaded() const;
	void resetOverloadCount()    { overloadCount = 0;      }
	int incrementOverloadCount() { return overloadCount++; }
	int getOverloadCount() const { return overloadCount;   }

	void updateBasisPointers(std::vector<Basisfunction*> &basis) ;

	virtual void read(std::istream &is);
	virtual void write(std::ostream &os) const;

private:
	std::vector<double> min;  // lower left corner in typical 2 or 3 dimensions
	std::vector<double> max;  // upper right corner 
	int id_;

	HashSet<Basisfunction*> support_;
	std::vector<int> support_ids_; // temporary storage for the read() method only

	int overloadCount ;
	
};

} // end namespace LR

#endif

