#include "LRSpline/Basisfunction.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Meshline.h"
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <set>

typedef unsigned int uint;

namespace LR {

/************************************************************************************************************************//**
 * \brief Default bivariate constructor
 * \param dim The dimension in the physical space, i.e. the number of components of the controlpoints
 * \param order_u Polynomial order (degree + 1) in first parametric direction
 * \param order_v Polynomial order (degree + 1) in second parametric direction
 ***************************************************************************************************************************/
Basisfunction::Basisfunction(int dim, int order_u, int order_v) {
	weight_       = 1;
	id_           = -1;
	hashCode_     = 0;
	knots_.resize(2);
	knots_[0].resize(order_u+1);
	knots_[1].resize(order_v+1);
	controlpoint_.resize(dim);
}

/************************************************************************************************************************//**
 * \brief Destructor. Cleans up back-pointers from elements that point to this B-spline
 ***************************************************************************************************************************/
Basisfunction::~Basisfunction() {
	PROFILE("Function destruction");
	for(uint i=0; i<support_.size(); i++)
		support_[i]->removeSupportFunction( (Basisfunction*) this);
}

#ifdef HAS_GOTOOLS
/************************************************************************************************************************//**
 * \brief Get parametric greville parameter for this B-spline
 * \returns The internal knot average
 ***************************************************************************************************************************/
Go::Point Basisfunction::getGrevilleParameter() const {
	Go::Point ans(knots_.size());
	for (size_t d = 0; d < knots_.size(); ++d) {
		ans[d] = 0.0;
		for(uint i=1; i<knots_[d].size()-1; i++)
			ans[d] += knots_[d][i];
		ans[d] /= (knots_[d].size()-2);
	}
	return ans;
}

/************************************************************************************************************************//**
 * \brief Get the control point
 * \param pt [out] The ascociated control point to this B-spline
 ***************************************************************************************************************************/
void Basisfunction::getControlPoint(Go::Point &pt) const {
	pt.resize(controlpoint_.size());
	for(uint d=0; d<controlpoint_.size(); d++)
		pt[d] = controlpoint_[d];
}
#endif


/************************************************************************************************************************//**
 * \brief Get parametric greville parameter for this B-spline
 * \returns The internal knot average
 ***************************************************************************************************************************/
void Basisfunction::getGrevilleParameter(std::vector<double> &pt) const {
	pt.resize(knots_.size());
	for(uint i=0; i<knots_.size(); i++) {
		pt[i] = 0;
		for(uint j=1; j<knots_[i].size()-1; j++)
			pt[i] += knots_[i][j];
		pt[i] /= (knots_[i].size()-2);
	}
}


#if 0
// this is left in as perhaps a more readable verison of what's going on in the evaluation stuff

void Basisfunction::evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right, bool v_from_right) const {
	results.resize((derivs+1)*(derivs+2)/2);
	fill(results.begin(), results.end(), 0);
	if(knots_[0][0] > u || u > knots_[0].back())
		return ;
	if(knots_[1][0] > v || v > knots_[1].back())
		return ;

	double ans_u[knots_[0].size()-1];
	double ans_v[knots_[1].size()-1];
	// double diff_u = 0;
	// double diff_v = 0;
	// double diff_2u[3];
	// double diff_2v[3];
	double **diff_u;
	double **diff_v;
	diff_u = new double*[derivs+1];
	diff_v = new double*[derivs+1];
	int diff_level;

	for(uint i=0; i<knots_[0].size()-1; i++) {
		if(u_from_right)
			ans_u[i] = (knots_[0][i] <= u && u <  knots_[0][i+1]) ? 1 : 0;
		else
			ans_u[i] = (knots_[0][i] <  u && u <= knots_[0][i+1]) ? 1 : 0;
	}

	diff_level = knots_[0].size()-2;
	for(uint n=1; n<knots_[0].size()-1; n++, diff_level--) {
		if(diff_level <= derivs) {
			diff_u[diff_level] = new double[diff_level+1];
			for(int j=0; j<=diff_level; j++)
				diff_u[diff_level][j] = ans_u[j];
		}
		for(int d = diff_level; d<= derivs; d++) {
			for(uint j=0; j<knots_[0].size()-1-n; j++) {
				diff_u[d][j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (   n   )/(knots_[0][j+n]  -knots_[0][ j ])*diff_u[d][ j ];
				diff_u[d][j] -= (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (   n   )/(knots_[0][j+n+1]-knots_[0][j+1])*diff_u[d][j+1];
			}
		}

#if 0
// kept this to maybe make it easier to see the logic behind evaluating arbitrary high derivatives
		if(n==knots_[0].size()-2)
			for(int j=0; j<=knots_[0].size()-n; j++)
				diff_2u[j] = ans_u[j];
		if(n>=knots_[0].size()-2) {
			for(int j=0; j<knots_[0].size()-n; j++) {
				diff_2u[j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (   n   )/(knots_[0][j+n]  -knots_[0][ j ])*diff_2u[ j ];
				diff_2u[j] -= (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (   n   )/(knots_[0][j+n+1]-knots_[0][j+1])*diff_2u[j+1];
			}
		}
		if(n==knots_[0].size()-1) {
			int j=0;
			diff_u  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (   knots_[0].size()-1   )/(knots_[0][j+n]  -knots_[0][ j ])*ans_u[ j ];
			diff_u -= (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (   knots_[0].size()-1   )/(knots_[0][j+n+1]-knots_[0][j+1])*ans_u[j+1];
		}
#endif
		for(uint j=0; j<knots_[0].size()-1-n; j++) {
			ans_u[j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (  u-knots_[0][j]  )/(knots_[0][j+n]  -knots_[0][ j ])*ans_u[ j ];
			ans_u[j] += (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (knots_[0][j+n+1]-u)/(knots_[0][j+n+1]-knots_[0][j+1])*ans_u[j+1];
		}
	}

	for(uint i=0; i<knots_[1].size()-1; i++) {
		if(v_from_right)
			ans_v[i] = (knots_[1][i] <= v && v <  knots_[1][i+1]) ? 1 : 0;
		else
			ans_v[i] = (knots_[1][i] <  v && v <= knots_[1][i+1]) ? 1 : 0;
	}

	diff_level = knots_[1].size()-2;
	for(uint n=1; n<knots_[1].size()-1; n++, diff_level--) {
		if(diff_level <= derivs) {
			diff_v[diff_level] = new double[diff_level+1];
			for(int j=0; j<=diff_level; j++)
				diff_v[diff_level][j] = ans_v[j];
		}
		for(int d = diff_level; d<= derivs; d++) {
			for(uint j=0; j<knots_[1].size()-1-n; j++) {
				diff_v[d][j]  = (knots_[1][ j+n ]==knots_[1][ j ]) ? 0 : (   n   )/(knots_[1][j+n]  -knots_[1][ j ])*diff_v[d][ j ];
				diff_v[d][j] -= (knots_[1][j+n+1]==knots_[1][j+1]) ? 0 : (   n   )/(knots_[1][j+n+1]-knots_[1][j+1])*diff_v[d][j+1];
			}
		}

		for(uint j=0; j<knots_[1].size()-1-n; j++) {
			ans_v[j]  = (knots_[1][ j+n ]==knots_[1][ j ]) ? 0 : (  v-knots_[1][j]  )/(knots_[1][j+n]  -knots_[1][ j ])*ans_v[ j ];
			ans_v[j] += (knots_[1][j+n+1]==knots_[1][j+1]) ? 0 : (knots_[1][j+n+1]-v)/(knots_[1][j+n+1]-knots_[1][j+1])*ans_v[j+1];
		}
	}

	// collect results
	diff_u[0] = ans_u;
	diff_v[0] = ans_v;
	int ip=0;
	for(int totDeriv=0; totDeriv<=derivs; totDeriv++)
		for(int vDeriv=0; vDeriv<=totDeriv; vDeriv++)
			results[ip++] = diff_u[totDeriv-vDeriv][0] * diff_v[vDeriv][0] * weight_;

	// clean up memory allocations
	for(int i=1; i<derivs+1; i++) {
		delete[] diff_u[i];
		delete[] diff_v[i];
	}
	delete[] diff_u;
	delete[] diff_v;

}

double Basisfunction::evaluate(double u, double v, bool u_from_right, bool v_from_right) const {
	if(knots_[0][0] > u || u > knots_[0].back())
		return 0;
	if(knots_[1][0] > v || v > knots_[1].back())
		return 0;

	double ans_u[knots_[0].size()-1];
	double ans_v[knots_[1].size()-1];

	for(uint i=0; i<knots_[0].size()-1; i++) {
		if(u_from_right)
			ans_u[i] = (knots_[0][i] <= u && u <  knots_[0][i+1]) ? 1 : 0;
		else
			ans_u[i] = (knots_[0][i] <  u && u <= knots_[0][i+1]) ? 1 : 0;
	}
	for(uint n=1; n<knots_[0].size()-1; n++)
		for(uint j=0; j<knots_[0].size()-1-n; j++) {
			ans_u[j]  = (knots_[0][ j+n ]==knots_[0][ j ]) ? 0 : (  u-knots_[0][j]  )/(knots_[0][j+n]  -knots_[0][ j ])*ans_u[ j ];
			ans_u[j] += (knots_[0][j+n+1]==knots_[0][j+1]) ? 0 : (knots_[0][j+n+1]-u)/(knots_[0][j+n+1]-knots_[0][j+1])*ans_u[j+1];
	}

	for(uint i=0; i<knots_[1].size()-1; i++) {
		if(v_from_right)
			ans_v[i] = (knots_[1][i] <= v && v <  knots_[1][i+1]) ? 1 : 0;
		else
			ans_v[i] = (knots_[1][i] <  v && v <= knots_[1][i+1]) ? 1 : 0;
	}
	for(uint n=1; n<knots_[1].size()-1; n++)
		for(uint j=0; j<knots_[1].size()-1-n; j++) {
			ans_v[j]  = (knots_[1][ j+n ]==knots_[1][ j ]) ? 0 : (  v-knots_[1][j]  )/(knots_[1][j+n]  -knots_[1][ j ])*ans_v[ j ];
			ans_v[j] += (knots_[1][j+n+1]==knots_[1][j+1]) ? 0 : (knots_[1][j+n+1]-v)/(knots_[1][j+n+1]-knots_[1][j+1])*ans_v[j+1];
	}

	return ans_u[0]*ans_v[0]*weight_;
}
#endif

/************************************************************************************************************************//**
 * \brief evaluates a bivariate B-spline
 * \param u Parametric evaluation point
 * \param v Parametric evaluation point
 * \param u_from_right Evaluate first parametric coordinate in the limit from the right
 * \param v_from_right Evaluate second parametric coordinate in the limit from the right
 * \return The B-spline evaluated at the chosen parametric coordinate
 ***************************************************************************************************************************/
double Basisfunction::evaluate(double u, double v, bool u_from_right, bool v_from_right) const {
	std::vector<double> results;
	std::vector<double> parPt(2);
	std::vector<bool>   fromRight(2);
	parPt[0] = u;
	parPt[1] = v;
	fromRight[0] = u_from_right;
	fromRight[1] = v_from_right;
	evaluate(results, parPt, 0, fromRight);
	return results[0];
}

/************************************************************************************************************************//**
 * \brief evaluates a trivariate B-spline
 * \param u Parametric evaluation point
 * \param v Parametric evaluation point
 * \param w Parametric evaluation point
 * \param u_from_right Evaluate first parametric coordinate in the limit from the right
 * \param v_from_right Evaluate second parametric coordinate in the limit from the right
 * \param w_from_right Evaluate third parametric coordinate in the limit from the right
 * \return The B-spline evaluated at the chosen parametric coordinate
 ***************************************************************************************************************************/
double Basisfunction::evaluate(double u, double v, double w, bool u_from_right, bool v_from_right, bool w_from_right) const {
	std::vector<double> results;
	std::vector<double> parPt(3);
	std::vector<bool>   fromRight(3);
	parPt[0] = u;
	parPt[1] = v;
	parPt[2] = w;
	fromRight[0] = u_from_right;
	fromRight[1] = v_from_right;
	fromRight[2] = w_from_right;
	evaluate(results, parPt, 0, fromRight);
	return results[0];
}

/************************************************************************************************************************//**
 * \brief evaluates a bivariate B-spline
 * \param results [out] Vector of all results. Upon function return contains (derivs+1)*(derivs+2)/2 evaluations of the B-spline
 *                      itself and all possible cross-derivatives up to order derivs. These are ordered as 1,dx,dy,d2x,dxdy,d2y,...
 * \param u Parametric evaluation point
 * \param v Parametric evaluation point
 * \param derivs Number of derivatives requested
 * \param u_from_right Evaluate first parametric coordinate in the limit from the right
 * \param v_from_right Evaluate second parametric coordinate in the limit from the right
 ***************************************************************************************************************************/
void Basisfunction::evaluate(std::vector<double> &results, double u, double v, int derivs, bool u_from_right, bool v_from_right) const {
	std::vector<double> parPt(2);
	std::vector<bool>   fromRight(2);
	parPt[0] = u;
	parPt[1] = v;
	fromRight[0] = u_from_right;
	fromRight[1] = v_from_right;
	evaluate(results, parPt, derivs, fromRight);
}

/************************************************************************************************************************//**
 * \brief evaluates a trivariate B-spline
 * \param results [out] Vector of all results. Upon function return contains (derivs+1)*(derivs+2)*(2*derivs+6)/12 evaluations of the B-spline
 *                      itself and all possible cross-derivatives up to order derivs. These are ordered as 1,dx,dy,d2x,dxdy,dxdz,d2y,...
 * \param u Parametric evaluation point
 * \param v Parametric evaluation point
 * \param w Parametric evaluation point
 * \param derivs Number of derivatives requested
 * \param u_from_right Evaluate first parametric coordinate in the limit from the right
 * \param v_from_right Evaluate second parametric coordinate in the limit from the right
 * \param w_from_right Evaluate third parametric coordinate in the limit from the right
 ***************************************************************************************************************************/
void Basisfunction::evaluate(std::vector<double> &results, double u, double v, double w, int derivs, bool u_from_right, bool v_from_right, bool w_from_right) const {
	std::vector<double> parPt(3);
	std::vector<bool>   fromRight(3);
	parPt[0] = u;
	parPt[1] = v;
	parPt[2] = w;
	fromRight[0] = u_from_right;
	fromRight[1] = v_from_right;
	fromRight[2] = w_from_right;
	evaluate(results, parPt, derivs, fromRight);
}

// funky recursive algorithm for collecting derivative results. Best illustrated by examples:
// ordering for bivariate second derivatives:  1, dx,dy, d2x,dxdy,d2y
// ordering for trivariate second derivatives: 1, dx,dy,dz, d2x,dxdy,dxdz,d2y,dydz,d2z
// ordering for trivariate third derivatives:  1,
//                                            dx,dy,dz,
//                                            d2x,dxdy,dxdz,d2y,dydz,d2z,
//                                            d3x,d2xdy,d2xdz,dxd2y,dxdydz,dxd2z,d3y,d2ydz,dyd2z,d3z
void collectResults(std::vector<double>::iterator &result,
                    double product,
                    std::vector<std::vector<std::vector<double> > > &diff,
                    int derivsLeft,
                    uint dim) {
	if(dim == diff.size()-1) {
		*result *= product*diff.back()[derivsLeft][0];
		result++;
		return;
	}
	for(int d=derivsLeft; d>-1; d--) {
		double ans = diff[dim][d][0];
		collectResults(result, ans*product, diff, derivsLeft-d, dim+1);
	}
};

/************************************************************************************************************************//**
 * \brief evaluates a general B-spline (currently only bivariate and trivariate supported - small fix to extend, but not now)
 * \param results [out] Vector of all results
 * \param parPt Parametric evaluation point
 * \param derivs Number of derivatives requested
 * \param from_right Vector of same size as parPt stating if any of the parametric directions should be evaluated in the limit from
 *                   the right
 *
 * Upon function return, results contains either (derivs+1)*(derivs+2)/2 (bivariate) or (derivs+1)*(derivs+2)*(2*derivs+6)/12 (trivariate)
 * evaluation points. These are all derivatives and organized as follows:
 * Bivariate splines up to second order: 1,dx,dy,d2x,dxdy,d2y,...
 * Trivariate splines up to second order: 1, dx,dy,dz, d2x,dxdy,dxdz,d2y,dydz,d2z
 * Trivariate splines up to third order: 1, dx,dy,dz, d2x,dxdy,dxdz,d2y,dydz,d2z, d3x,d2xdy,d2xdz,dxd2y,dxdydz,dxd2z,d3y,d2ydz,dyd2z,d3z
 ***************************************************************************************************************************/
void Basisfunction::evaluate(std::vector<double> &results, const std::vector<double> &parPt, int derivs, const std::vector<bool> &from_right) const {
	uint dim = knots_.size();
	if(dim != parPt.size() || dim != from_right.size()) {
		std::cerr << "Error Basisfunction::evalate(...) parametric dimension mismatch" << std::endl;
		exit(9230);
	}

	if(dim == 2) {    // bivariate splines
		results.resize((derivs+1)*(derivs+2)/2);               // (this is the triangular numbers)
	} else if(dim == 3) { // trivariate splines
		results.resize((derivs+1)*(derivs+2)*(2*derivs+6)/12); // (sum of triangular numbers)
	} else {
		std::cerr << "Error Basisfunction::evalate(...) for parametric dimension other than 2 or 3" << std::endl;
		exit(9231);
	}
	fill(results.begin(), results.end(), 0.0);

	std::vector<std::vector<double> >               ans(dim);
	std::vector<std::vector<std::vector<double> > > diff(dim);
	uint i = 0;
	for(std::vector<double> knot : knots_) {
		if(knot[0] > parPt[i] || parPt[i] > knot.back())
			return;
		ans[i].resize(knot.size()-1);
		for(uint j=0; j<knot.size()-1; j++) {
			if(from_right[i])
				ans[i][j] = (knot[j] <= parPt[i] && parPt[i] <  knot[j+1]) ? 1 : 0;
			else
				ans[i][j] = (knot[j] <  parPt[i] && parPt[i] <= knot[j+1]) ? 1 : 0;
		}

		int p          = knot.size()-2;
		int diff_level = p;
		diff[i].resize(derivs+1, std::vector<double>(1, 0));
		for(uint n=1; n<knot.size()-1; n++, diff_level--) {
			if(diff_level <= derivs) {
				diff[i][diff_level].resize(diff_level+1);
				for(int j=0; j<=diff_level; j++)
					diff[i][diff_level][j] = ans[i][j];
			}
			for(int d = diff_level; d <= derivs && d <= p; d++) {
				for(uint j=0; j<knot.size()-1-n; j++) {
					diff[i][d][j]  = (knot[ j+n ]==knot[ j ]) ? 0 : (   n   )/(knot[j+n]  -knot[ j ])*diff[i][d][ j ];
					diff[i][d][j] -= (knot[j+n+1]==knot[j+1]) ? 0 : (   n   )/(knot[j+n+1]-knot[j+1])*diff[i][d][j+1];
				}
			}
			for(uint j=0; j<knot.size()-1-n; j++) {
				ans[i][j]  = (knot[ j+n ]==knot[ j ]) ? 0 : (  parPt[i]-knot[j]  )/(knot[j+n]  -knot[ j ])*ans[i][ j ];
				ans[i][j] += (knot[j+n+1]==knot[j+1]) ? 0 : (knot[j+n+1]-parPt[i])/(knot[j+n+1]-knot[j+1])*ans[i][j+1];
			}
		}

		i++;
	}

	// collect results
	for(i=0; i<dim; i++)
		diff[i][0] = ans[i];

	fill(results.begin(), results.end(), weight_);
	std::vector<double>::iterator resIt = results.begin();
	for(int totDeriv=0; totDeriv<=derivs; totDeriv++)
		collectResults(resIt, 1.0, diff, totDeriv, 0);
}

/************************************************************************************************************************//**
 * \brief Get the control point
 * \param pt [out] The ascociated control point to this B-spline
 ***************************************************************************************************************************/
void Basisfunction::getControlPoint(std::vector<double> &pt) const {
	pt.resize(controlpoint_.size());
	for(uint d=0; d<controlpoint_.size(); d++)
		pt[d] = controlpoint_[d];
}

/************************************************************************************************************************//**
 * \brief Remove element from the list of supported elements
 * \param el The element to remove
 * \return True if the element was found and successfully removed
 ***************************************************************************************************************************/
bool Basisfunction::removeSupport(Element *el) {
	for(uint i=0; i<support_.size(); i++) {
		if(el == support_[i]) {
			support_[i] = support_.back();
			support_.pop_back();
			return true;
		}
	}
	return false;
}

/************************************************************************************************************************//**
 * \brief Add element from the list of supported elements if its support overlaps
 * \param el The element to add
 * \return True if the elements support overlaps with the support of this B-spline
 ***************************************************************************************************************************/
bool Basisfunction::addSupport(Element *el) {
	if(overlaps(el)) {
		support_.push_back(el);
		return true;
	}
	return false;
}

/************************************************************************************************************************//**
 * \brief Returns true if this B-splines support overlaps with the elements size
 ***************************************************************************************************************************/
bool Basisfunction::overlaps(Element *el) const {
	for(uint i=0; i<knots_.size(); i++) {
		if(knots_[i][0]     >= el->getParmax(i))
			return false;
		if(knots_[i].back() <= el->getParmin(i))
			return false;
	}
	return true;
}

/************************************************************************************************************************//**
 * \brief Get all B-splines which have overlapping support with this B-spline. The combined support of this set is denoted
 *        the extended support of the B-spline
 * \return A list of elements which describes the extended support
 ***************************************************************************************************************************/
std::vector<Element*> Basisfunction::getExtendedSupport() {
	std::set<Element*> ans ;

	for(Element* e : support_)
		for(Basisfunction* b : e->support())
			for(Element* e2 : b->support())
				ans.insert(e2);
	std::vector<Element*> ans_vector(ans.size());
	std::copy(ans.begin(), ans.end(), ans_vector.begin());
	return ans_vector;
}

/************************************************************************************************************************//**
 * \brief Get a larger support for the B-splines. Typically trying our best to just extend the support by one element in one direction
 * \return A list of elements which describes the minimal extended support
 ***************************************************************************************************************************/
std::vector<Element*> Basisfunction::getMinimalExtendedSupport() {
	if(knots_.size() != 2) {
		std::cerr << "Error: Basisfunction::getMinimalExtendedSupport() only for bivariate B-splines" << std::endl;
		exit(86136);
	}
	double min_du = DBL_MAX;
	double min_dv = DBL_MAX;
	Basisfunction *smallestGuy = NULL;

	bool edgeUmin = (knots_[0][0] == knots_[0][knots_[0].size()-2]);
	bool edgeUmax = (knots_[0][1] == knots_[0][knots_[0].size()-1]);
	bool edgeVmin = (knots_[1][0] == knots_[1][knots_[1].size()-2]);
	bool edgeVmax = (knots_[1][1] == knots_[1][knots_[1].size()-1]);

	if(! (edgeUmin || edgeUmax) )
		min_du = getParmax(0) - getParmin(0);
	if(! (edgeVmin || edgeVmax) )
		min_dv = getParmax(1) - getParmin(1);

	for(Element* e : support_) {
		for(Basisfunction* b : e->support()) {
			if( ((b->getParmin(0) <  this->getParmin(0) && b->getParmax(0) >= this->getParmax(0)) ||
			     (b->getParmin(0) <= this->getParmin(0) && b->getParmax(0) >  this->getParmax(0))) &&  // extend left OR right of the u-directions
			    ((b->getParmin(1) <  this->getParmin(1) && b->getParmax(1) >= this->getParmax(1)) ||
			     (b->getParmin(1) <= this->getParmin(1) && b->getParmax(1) >  this->getParmax(1))) &&  // AND extend up OR down in the v-direction
			    min_du >= b->getParmax(0)-b->getParmin(0) &&
			    min_dv >= b->getParmax(1)-b->getParmin(1)  ) { // AND less support than the last best guy

				min_du = b->getParmax(0)-b->getParmin(0);
				min_dv = b->getParmax(1)-b->getParmin(1);
				smallestGuy = b;
			}
		}
	}
	if(smallestGuy == NULL) {
		for(Element* e : support_) {
			for(Basisfunction* b : e->support()) {
				if(  b->getParmin(0) <= this->getParmin(0)    &&
				     b->getParmax(0) >= this->getParmax(0)    &&
				     b->getParmin(1) <= this->getParmin(1)    &&
				     b->getParmax(1) >= this->getParmax(1)    &&  // Cover *this' support
				    (b->getParmin(0) <  this->getParmin(0) ||
				     b->getParmax(0) >  this->getParmax(0) ||
				     b->getParmin(1) <  this->getParmin(1) ||    // Extend at least one direction
				     b->getParmax(1) >  this->getParmax(1) )  &&
					min_du >= b->getParmax(0)-b->getParmin(0) &&
					min_dv >= b->getParmax(1)-b->getParmin(1)  ) { // smaller support than the last best option

					min_du = b->getParmax(0)-b->getParmin(0);
					min_dv = b->getParmax(1)-b->getParmin(1);
					smallestGuy = b;
				}
			}
		}
	}

	// Not all B-splines do have a proper minimal extended support as defined above
	if(smallestGuy == NULL)
		return std::vector<Element*>(0);

	std::vector<Element*> results(smallestGuy->nSupportedElements());
	std::copy(smallestGuy->supportedElementBegin(), smallestGuy->supportedElementEnd(), results.begin());
	return results;
}


/************************************************************************************************************************//**
 * \brief Returns all functions which has support on one or more of the elements that this function have support on
 ***************************************************************************************************************************/
HashSet<Basisfunction*> Basisfunction::getOverlappingFunctions() const {
	HashSet<Basisfunction*> result;
	for(auto e : support())        // for all elements that this have support on
		for(auto b : e->support()) // for all basis functions with support on that element
			result.insert(b);
	return result;
}



/************************************************************************************************************************//**
 * \brief Resize the B-spline dimension. All control point values set to 0.0
 * \param dim New control point dimension
 ***************************************************************************************************************************/
void Basisfunction::setDimension(int dim) {
	controlpoint_.resize(dim);
	for(int i=0; i<dim; i++)
		controlpoint_[i]  = 0.0;
}

/************************************************************************************************************************//**
 * \brief Get the B-spline hash code for storage in the HashSet container
 * \returns some "random" long based on the local knot vector
 ***************************************************************************************************************************/
long Basisfunction::hashCode() const {
	if(hashCode_ != 0)
		return hashCode_;

	int nKnots = 0;
	for(std::vector<double> knot : knots_)
		nKnots += knot.size();

	int bitsFromEach = (sizeof(long)*8) / nKnots;
	int bitsLeft     = (sizeof(long)*8) % nKnots;
	int offset       = 0;
	for(std::vector<double> knot : knots_) {
		for(uint i=0; i<knot.size()-1; i++) {
			long randInt = log2(fabs(knot[i]))*120000;
			// int randInt = log2(fabs(knot[i]+1))*343;
			// int randInt = (i%2==0) ? knot[i] : log2(fabs(knot[i]));
			if(knot[i]==0) randInt = 0;
			// int randInt = log2(fabs(knot[i]));
			// int randInt = knot[i];
			long mask    = 0;
			for(int k=0; k<bitsFromEach + (bitsLeft>0); k++)
				mask |= (1<<k);
			hashCode_ |= ((long) (randInt & mask)) << offset;
			offset += bitsFromEach + (bitsLeft>0);
			bitsLeft--;
		}
	}
	return hashCode_;
}

/************************************************************************************************************************//**
 * \brief Test for B-spline equality
 * \param other The other B-spline to check against
 * \returns True if the knot vectors are identical (up to a tolerance of 1e-10)
 ***************************************************************************************************************************/
bool Basisfunction::equals(const Basisfunction &other) const {
	if(knots_.size() != other.knots_.size())
		return false;
	for(uint i=0; i<knots_.size(); i++) {
		if(knots_[i].size() != other[i].size())
			return false;
		for(uint j=0; j<knots_[i].size(); j++)
			if(fabs(knots_[i][j] - other[i][j]) > 1e-10)
				return false;
	}
	return true;
}

/************************************************************************************************************************//**
 * \brief Test for B-spline equality
 * \param other The other B-spline to check against
 * \returns True if the knot vectors are identical (up to a tolerance of 1e-10)
 ***************************************************************************************************************************/
bool Basisfunction::operator==(const Basisfunction &other) const {
	return equals(other);
}

/************************************************************************************************************************//**
 * \brief Test for B-spline equality
 * \param other The other B-spline to check against
 * \returns True if the knot vectors are identical (up to a tolerance of 1e-10)
 ***************************************************************************************************************************/
void Basisfunction::operator+=(const Basisfunction &other) {
	double newWeight = weight_ + other.weight_;
	for(uint i=0; i<controlpoint_.size(); i++)
		controlpoint_[i] = (controlpoint_[i]*weight_ + other.controlpoint_[i]*other.weight_)/newWeight;
	weight_ = newWeight;
}

/************************************************************************************************************************//**
 * \brief Returns a deep copy of the B-spline with all the same knot vectors and controlpoints. Note: Does not copy supported elements
 ***************************************************************************************************************************/
Basisfunction* Basisfunction::copy() const {

	std::vector<int> order;
	for(uint i=0; i<knots_.size(); i++)
		order.push_back(knots_[i].size()-1);
	Basisfunction *returnValue = new Basisfunction(controlpoint_.size(), knots_.size(), order);

	for(uint i=0; i<knots_.size(); i++)
		std::copy(knots_[i].begin(), knots_[i].end(), returnValue->knots_[i].begin());

	std::copy(controlpoint_.begin(), controlpoint_.end(), returnValue->controlpoint_.begin());
	returnValue->weight_ = weight_;
	returnValue->id_     = id_;

	return returnValue;
}




/************************************************************************************************************************//**
 * \brief Reads a B-spline from input stream
 ***************************************************************************************************************************/
void Basisfunction::read(std::istream &is) {
// convenience macro for reading formated input
#define ASSERT_NEXT_CHAR(c) {ws(is); nextChar = is.get(); if(nextChar!=c) { std::cerr << "Error parsing basis function\n"; exit(324); } ws(is); }
	char nextChar;

	// read id tag
	is >> id_;
	ws(is);
	ASSERT_NEXT_CHAR(':');

	// read knot vectors
	bool isFirst = true;
	for(uint i=0; i<knots_.size(); i++) {
		if(!isFirst) ASSERT_NEXT_CHAR('x');
		ASSERT_NEXT_CHAR('[');
		for(uint j=0; j<knots_[i].size(); j++)
			is >> knots_[i][j];
		ASSERT_NEXT_CHAR(']');
		isFirst = false;
	}

	// read control point
	for(uint i=0; i<controlpoint_.size(); i++)
		is >> controlpoint_[i];

	// read weight
	ASSERT_NEXT_CHAR('(');
	is >> weight_;
	ASSERT_NEXT_CHAR(')');
#undef ASSERT_NEXT_CHAR
}

/************************************************************************************************************************//**
 * \brief Writes a B-spline to output stream
 ***************************************************************************************************************************/
void Basisfunction::write(std::ostream &os) const {
	os << id_ << ": ";
	bool isFirst = true;
	for(std::vector<double> knot : knots_) {
		if(!isFirst) os << "x ";
		os << "[";
		for(uint i=0; i<knot.size(); i++)
			os << knot[i] << " ";
		os << "] ";
		isFirst = false;
	}

	for(uint i=0; i<controlpoint_.size(); i++)
		os << controlpoint_[i] << " ";
	os << "(" << weight_ << ")";
}

/************************************************************************************************************************//**
 * \brief Returns true if all supported elements are considered overloaded (i.e. has more than (p+1)*(p+1) supported B-splines)
 ***************************************************************************************************************************/
bool Basisfunction::isOverloaded()  const {
	for(uint i=0; i<support_.size(); i++)
		if(! support_[i]->isOverloaded() )
			return false;
	return true;
}

/************************************************************************************************************************//**
 * \brief Returns the maximum overloaded count of the supported elements
 ***************************************************************************************************************************/
int Basisfunction::getOverloadCount() const {
	int ans = INT_MAX;
	for(uint i=0; i<support_.size(); i++)
		ans = (ans < support_[i]->getOverloadCount()) ? ans : support_[i]->getOverloadCount();
	return ans;
}

/************************************************************************************************************************//**
 * \brief Returns true if the support of *other is completely contained in *this
 ***************************************************************************************************************************/
bool Basisfunction::contains(const Basisfunction &other) const {
	for(uint i=0; i<knots_.size(); i++)
		if(other[i][0]     < (*this)[i][0] ||
		   other[i].back() > (*this)[i].back())
			return false;
	return true;
}

/************************************************************************************************************************//**
 * \brief flip two parametric coordinate directions
 ***************************************************************************************************************************/
void Basisfunction::flip(int dir1, int dir2) {
	std::vector<double> tmp = knots_[dir1];
	knots_[dir1] = knots_[dir2];
	knots_[dir2] = tmp;
}

/************************************************************************************************************************//**
 * \brief reverse one parametric direction. Need global range (parmin, parmax) for scaling
 ***************************************************************************************************************************/
void Basisfunction::reverse(int pardir, double parmin, double parmax) {
	std::vector<double> tmp(knots_[pardir]);
	int n = tmp.size();
	for(int i=0; i<n; i++) {
		knots_[pardir][n-i-1] = (parmax - tmp[i]) / (parmax-parmin) * (parmax-parmin) + parmin;
	}
}

/************************************************************************************************************************//**
 * \brief (used when iterating over all functions), scales knot vectors so they globally fit into range (0,1)
 ***************************************************************************************************************************/
void Basisfunction::normalize(int pardir, double parmin, double parmax) {
	for(uint i=0; i<knots_[pardir].size(); i++) {
		knots_[pardir][i] = (knots_[pardir][i] - parmin) / (parmax-parmin) * (parmax-parmin) + parmin;
	}
}


} // end namespace LR
