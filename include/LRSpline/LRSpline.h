#ifndef LRSPLINE_H
#define LRSPLINE_H

#include "HashSet.h"
#include "HashSet.h"
#include <GoTools/geometry/Streamable.h>
#include <vector>

enum refinementStrategy {
	LR_SAFE = 0,
	LR_MINSPAN = 1,
	LR_ISOTROPIC_EL = 2,
	LR_ISOTROPIC_FUNC = 3,
};


namespace LR {

// forms a 6 bit binary mask where each bit corresponds to a side
// corners lines have 2 bits 'on'
// corners vertices have 3 bits 'on'
enum parameterEdge {
NONE       = 0,
WEST       = 1,    // 000001
EAST       = 2,    // 000010
SOUTH      = 4,    // 000100
NORTH      = 8,    // 001000
TOP        = 16,   // 010000
BOTTOM     = 32,   // 100000
// convienience variables for 2D case follows. In general use SOUTH|WEST, SOUTH|EAST, etc
SOUTH_WEST = 5,    // 000101
SOUTH_EAST = 7,    // 000110
NORTH_WEST = 9,    // 001001
NORTH_EAST = 10};  // 001010

typedef enum parameterEdge parameterEdge;

class Basisfunction;
class Element;

class LRSpline : public Go::Streamable {

public:
	LRSpline();

	virtual void generateIDs() const;

	// common get methods
	int nBasisFunctions()    const { return basis_.size()  ; };
	int nElements()          const { return element_.size(); };
	int dimension()          const { return dim_           ; };
	int order        (int i) const { return order_[i]      ; };
	double startparam(int i) const { return start_[i]      ; };
	double endparam  (int i) const { return end_[i]        ; };
	bool rational()          const { return rational_      ; };
	
	// more funky get methods
	void getEdgeFunctions(std::vector<Basisfunction*> &edgeFunctions, parameterEdge edge, int depth=1) const;
	void getEdgeElements( std::vector<Element*>       &edgeElements,  parameterEdge edge             ) const;

	// get container iterators
	std::vector<Element*>::iterator        elementBegin()         { return element_.begin(); };
	std::vector<Element*>::iterator        elementEnd()           { return element_.end();   };
	HashSet_iterator<Basisfunction*>       basisBegin()           { return basis_.begin();   };
	HashSet_iterator<Basisfunction*>       basisEnd()             { return basis_.end();     };
	HashSet_const_iterator<Basisfunction*> basisBegin()     const { return basis_.begin();   };
	HashSet_const_iterator<Basisfunction*> basisEnd()       const { return basis_.end();     };
	const HashSet<Basisfunction*>& getAllBasisfunctions()   const { return basis_ ;          };
	const std::vector<Element*>&           getAllElements() const { return element_ ;        };

	// refinement functions
	virtual void refineBasisFunction(int index) = 0;
	virtual void refineBasisFunction(std::vector<int> &indices) = 0;
	virtual void refineElement(int index) = 0;
	virtual void refineElement(const std::vector<int> &indices) = 0;
	virtual void refineByDimensionIncrease(const std::vector<double> &error, double beta) = 0;

	// set refinement state parameters
	void setRefStrat(enum refinementStrategy strat) { refStrat_        = strat;    };
	void setRefSymmetry(int symmetry)               { this->symmetry_  = symmetry; };
	void setRefMultiplicity(int mult)               { refKnotlineMult_ = mult;     };
	void setMaxTjoints(int n)                       { maxTjoints_      = n;        };
	void setCloseGaps(bool doClose)                 { doCloseGaps_     = doClose;  };
	void setMaxAspectRatio(double r, bool aposterioriFix=true) {
		maxAspectRatio_ = r;
		doAspectRatioFix_ = aposterioriFix;
	}

	// input output methods
	virtual void read(std::istream &is)         { };
	virtual void write(std::ostream &os) const  { };

protected:
	// useful descriptive stuff
	int dim_;
	bool rational_;
	std::vector<double> start_ ; //! \brief parametric start coordinate (2 components for surfaces, 3 for volumes)
	std::vector<double> end_   ; //! \brief parametric stop coordinate (2 components for surfaces, 3 for volumes)
	std::vector<int>    order_ ; //! \brief polynomial order (degree + 1) in each parametric direction (2 or 3 components)
	
	// core storage places for the building blocks
	std::vector<Basisfunction*> basisVector; // only used in read/write functions
	HashSet<Basisfunction*> basis_;
	std::vector<Element*> element_;

	// refinement parameters
	enum refinementStrategy refStrat_;
	int                     refKnotlineMult_;
	int                     symmetry_;
	int                     maxTjoints_;
	bool                    doCloseGaps_;
	bool                    doAspectRatioFix_;
	double                  maxAspectRatio_;
	
};

} // end namespace LR



#endif
