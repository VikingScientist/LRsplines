#ifndef LRSPLINE_H
#define LRSPLINE_H

#include "Basisfunction.h"
#include "HashSet.h"
#include <GoTools/geometry/Streamable.h>

enum refinementStrategy {
	LR_SAFE = 0,
	LR_MINSPAN = 1,
	LR_ISOTROPIC_EL = 2,
	LR_ISOTROPIC_FUNC = 3,
};

namespace LR {

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
