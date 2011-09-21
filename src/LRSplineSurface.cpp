
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Basisfunction.h"
#include "LRSpline/Meshline.h"
#include "LRSpline/Element.h"
#include "LRSpline/Profiler.h"
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <boost/rational.hpp>

typedef unsigned int uint;

namespace LR {


LRSplineSurface::LRSplineSurface() {
	rational_ = false;
	dim_      = 0;
	order_u_  = 0;
	order_v_  = 0;
	start_u_  = 0;
	start_v_  = 0;
	end_u_    = 0;
	end_v_    = 0;
	basis_    = std::vector<Basisfunction*>(0);
	meshline_ = std::vector<Meshline*>(0);
	element_  = std::vector<Element*>(0);
}

LRSplineSurface::LRSplineSurface(Go::SplineSurface *surf) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	rational_ = surf->rational();
	dim_      = surf->dimension();
	order_u_  = surf->order_u();
	order_v_  = surf->order_v();
	start_u_  = surf->startparam_u();
	start_v_  = surf->startparam_v();
	end_u_    = surf->endparam_u();
	end_v_    = surf->endparam_v();

	int n1 = surf->numCoefs_u();
	int n2 = surf->numCoefs_v();
	std::vector<double>::iterator coef  = surf->coefs_begin();
	std::vector<double>::const_iterator knot_u = surf->basis(0).begin();
	std::vector<double>::const_iterator knot_v = surf->basis(1).begin();

	basis_ = std::vector<Basisfunction*>(n1*n2);
	int k=0;
	for(int j=0; j<n2; j++) {
		for(int i=0; i<n1; i++) {
			basis_[k++] = new Basisfunction(&(*(knot_u+i)), &(*(knot_v+j)), &(*(coef+(j*n1+i)*(dim_+rational_))), dim_, order_u_, order_v_);
			if(i == 0)    basis_[k-1]->addEdge(WEST);
			if(i == n1-1) basis_[k-1]->addEdge(EAST);
			if(j == 0)    basis_[k-1]->addEdge(SOUTH);
			if(j == n2-1) basis_[k-1]->addEdge(NORTH);
		}
	}
	int unique_u=0;
	int unique_v=0;
	for(int i=0; i<n1+order_u_; i++) {// const u, spanning v
		int mult = 1;
		while(i<n1+order_u_ && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2], mult) );
	}
	for(int i=0; i<n2+order_v_; i++) {// const v, spanning u
		int mult = 1;
		while(i<n2+order_v_ && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n2], mult) );
	}
	for(int j=0; j<unique_v-1; j++) {
		for(int i=0; i<unique_u-1; i++) {
			double umin = meshline_[i]->const_par_;
			double vmin = meshline_[unique_u + j]->const_par_;
			double umax = meshline_[i+1]->const_par_;
			double vmax = meshline_[unique_u + j+1]->const_par_;
			element_.push_back(new Element(umin, vmin, umax, vmax));
		}
	}
	for(uint i=0; i<basis_.size(); i++)
		updateSupport(basis_[i]);
}

LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational) {
#ifdef TIME_LRSPLINE
	PROFILE("Constructor");
#endif
	rational_ = rational;
	dim_      = dim;
	order_u_  = order_u;
	order_v_  = order_v;
	start_u_  = knot_u[0];
	start_v_  = knot_v[0];
	end_u_    = knot_u[n1];
	end_v_    = knot_v[n2];

	basis_ = std::vector<Basisfunction*>(n1*n2);
	int k=0;
	for(int j=0; j<n2; j++) {
		for(int i=0; i<n1; i++) {
			basis_[k++] = new Basisfunction(knot_u+i, knot_v+j, coef+(j*n1+i)*(dim+rational), dim, order_u, order_v);
			if(i == 0)    basis_[k-1]->addEdge(WEST);
			if(i == n1-1) basis_[k-1]->addEdge(EAST);
			if(j == 0)    basis_[k-1]->addEdge(SOUTH);
			if(j == n2-1) basis_[k-1]->addEdge(NORTH);
		}
	}
	int unique_u=0;
	int unique_v=0;
	for(int i=0; i<n1+order_u; i++) {// const u, spanning v
		int mult = 1;
		while(i<n1+order_u && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		unique_u++;
		meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2], mult) );
	}
	for(int i=0; i<n2+order_v; i++) {// const v, spanning u
		int mult = 1;
		while(i<n2+order_v && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		unique_v++;
		meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n2], mult) );
	}
	for(int j=0; j<unique_v-1; j++) {
		for(int i=0; i<unique_u-1; i++) {
			double umin = meshline_[i]->const_par_;
			double vmin = meshline_[unique_u + j]->const_par_;
			double umax = meshline_[i+1]->const_par_;
			double vmax = meshline_[unique_u + j+1]->const_par_;
			element_.push_back(new Element(umin, vmin, umax, vmax));
		}
	}
	for(uint i=0; i<basis_.size(); i++)
		updateSupport(basis_[i]);
}

void LRSplineSurface::point(Go::Point &pt, double u, double v) const {
	Go::Point cp;
	double basis_ev;
	pt.resize(dim_);
	pt *= 0;
	for(uint i=0; i<basis_.size(); i++) {
		basis_[i]->getControlPoint(cp);
		basis_ev = basis_[i]->evaluate(u,v, u!=end_u_, v!=end_v_);
		pt += basis_ev*cp;
	}
}

void LRSplineSurface::point(std::vector<Go::Point> &pts, double u, double v, int derivs) const {
#ifdef TIME_LRSPLINE
	PROFILE("Point()");
#endif
	Go::Point cp;
	std::vector<double> basis_ev;

	// clear and resize output array (optimization may consider this an outside task)
	pts.resize((derivs+1)*(derivs+2)/2);
	for(uint i=0; i<pts.size(); i++) {
		pts[i].resize(dim_);
		pts[i] *= 0;
	}

	for(uint i=0; i<basis_.size(); i++) {
		basis_[i]->getControlPoint(cp);
		basis_[i]->evaluate(basis_ev, u,v, derivs, u!=end_u_, v!=end_v_);
		for(uint j=0; j<pts.size(); j++)
			pts[j] += basis_ev[j]*cp;
	}
}

void LRSplineSurface::computeBasis (double param_u, double param_v, Go::BasisDerivsSf2 & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	std::vector<double> values;
	std::vector<Basisfunction*>::const_iterator it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->supportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->supportEnd();
	result.prepareDerivs(param_u, param_v, 0, -1, (itStop-itStart));
	int i=0;
	for(it=itStart; it!=itStop; it++, i++) {
		(*it)->evaluate(values, param_u, param_v, 2, param_u!=end_u_, param_v!=end_v_);
		result.basisValues[i]    = values[0];
		result.basisDerivs_u[i]  = values[1];
		result.basisDerivs_v[i]  = values[2];
		result.basisDerivs_uu[i] = values[3];
		result.basisDerivs_uv[i] = values[4];
		result.basisDerivs_vv[i] = values[5];
	}
}

void LRSplineSurface::computeBasis (double param_u, double param_v, Go::BasisDerivsSf & result, int iEl ) const {
#ifdef TIME_LRSPLINE
	PROFILE("computeBasis()");
#endif
	std::vector<double> values;
	std::vector<Basisfunction*>::const_iterator it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->supportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->supportEnd();
	result.prepareDerivs(param_u, param_v, 0, -1, (itStop-itStart));
	int i=0;
	for(it=itStart; it!=itStop; it++, i++) {
		(*it)->evaluate(values, param_u, param_v, 1, param_u!=end_u_, param_v!=end_v_);
		result.basisValues[i]   = values[0];
		result.basisDerivs_u[i] = values[1];
		result.basisDerivs_v[i] = values[2];
	}
}


void LRSplineSurface::computeBasis(double param_u, double param_v, Go::BasisPtsSf & result, int iEl ) const {
	std::vector<Basisfunction*>::const_iterator it, itStop, itStart;
	itStart = (iEl<0) ? basis_.begin() : element_[iEl]->supportBegin();
	itStop  = (iEl<0) ? basis_.end()   : element_[iEl]->supportEnd();
	result.preparePts(param_u, param_v, 0, -1, (itStop-itStart));
	int i=0;
	for(it=itStart; it!=itStop; it++, i++)
		result.basisValues[i] = (*it)->evaluate(param_u, param_v, param_u!=end_u_, param_v!=end_v_);
}

int LRSplineSurface::getElementContaining(double u, double v) const {
	for(uint i=0; i<element_.size(); i++)
		if(element_[i]->umin() <= u && element_[i]->vmin() <= v) 
			if((u < element_[i]->umax() || (u == end_u_ && u <= element_[i]->umax())) && 
			   (v < element_[i]->vmax() || (v == end_v_ && v <= element_[i]->vmax())))
				return i;
		
	return -1;
}

void LRSplineSurface::refineElement(int index) {
	std::vector<int> ref_index(1, index);
	refineElement(ref_index);
}

void LRSplineSurface::refineElement(std::vector<int> index) {
	// span-u lines
	std::vector<double> start_u(index.size());
	std::vector<double> stop_u (index.size());
	std::vector<double> mid_v  (index.size());

	// span-v lines
	std::vector<double> start_v(index.size());
	std::vector<double> stop_v (index.size());
	std::vector<double> mid_u  (index.size());

	std::vector<Basisfunction*>::iterator it;
	for(uint i=0; i<index.size(); i++) {
		double umin = element_[index[i]]->umin();
		double umax = element_[index[i]]->umax();
		double vmin = element_[index[i]]->vmin();
		double vmax = element_[index[i]]->vmax();
		for(it=element_[index[i]]->supportBegin(); it<element_[index[i]]->supportEnd(); it++) {
			umin = (umin > (**it).umin()) ? (**it).umin() : umin;
			umax = (umax < (**it).umax()) ? (**it).umax() : umax;
			vmin = (vmin > (**it).vmin()) ? (**it).vmin() : vmin;
			vmax = (vmax < (**it).vmax()) ? (**it).vmax() : vmax;
		}
		start_u[i] = umin;
		stop_u[i]  = umax;
		mid_v[i]   = (element_[index[i]]->vmin() + element_[index[i]]->vmax())/2.0;

		start_v[i] = vmin;
		stop_v[i]  = vmax;
		mid_u[i]   = (element_[index[i]]->umin() + element_[index[i]]->umax())/2.0;
	}
	for(uint i=0; i<index.size(); i++) {
		this->insert_const_u_edge(mid_u[i], start_v[i], stop_v[i], 1);
		this->insert_const_v_edge(mid_v[i], start_u[i], stop_u[i], 1);
	}
}

void LRSplineSurface::insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity) {
	insert_line(true, u, start_v, stop_v, multiplicity);
}

void LRSplineSurface::insert_line(bool const_u, double const_par, double start, double stop, int multiplicity) {
	Meshline *newline;
#ifdef TIME_LRSPLINE
	PROFILE("insert_line()");
#endif
	{ // check if the line is an extension or a merging of existing lines
#ifdef TIME_LRSPLINE
	PROFILE("line verification");
#endif
	newline = new Meshline(!const_u, const_par, start, stop, multiplicity);
	if(multiplicity != 1) {
		std::cerr << "ERROR: LRSplineSurface::insert_const_u_edge() not supported for multiplicity != 1\n";
		std::cerr << "       requested line: " << *newline << std::endl;
	}
	for(uint i=0; i<meshline_.size(); i++) {
		if( meshline_[i]->is_spanning_u() != const_u && meshline_[i]->const_par_ == const_par && 
		    meshline_[i]->start_ <= stop && meshline_[i]->stop_ >= start)  {
			// if newline overlaps any existing ones (may be multiple existing ones)
			// let newline be the entire length of all merged and delete the unused ones
		    	if(meshline_[i]->start_ < start) newline->start_ = meshline_[i]->start_;
		    	if(meshline_[i]->stop_  > stop ) newline->stop_  = meshline_[i]->stop_;
			meshline_.erase(meshline_.begin() + i);
			i--;
		}
	}
	}

	int nOldFunctions     = basis_.size();  // number of functions before any new ones are inserted
	int nRemovedFunctions = 0;
	int nNewFunctions     = 0; // number of newly created functions eligable for testing in STEP 2

	{ // STEP 1: test EVERY function against the NEW meshline
#ifdef TIME_LRSPLINE
	PROFILE("STEP 1");
#endif
	{
#ifdef TIME_LRSPLINE
	PROFILE("S1-basissplit");
#endif
	for(int i=0; i<nOldFunctions-nRemovedFunctions; i++) {
		if(newline->splits(basis_[i]) && !newline->containedIn(basis_[i])) {
			nNewFunctions += split( const_u, i, const_par );
			i--; // splitting deletes a basisfunction in the middle of the basis_ vector
			nRemovedFunctions++;
		}
	}
	}
	{
#ifdef TIME_LRSPLINE
	PROFILE("S1-elementsplit");
#endif
	for(uint i=0; i<element_.size(); i++) {
		if(newline->splits(element_[i]))
			element_.push_back(element_[i]->split(!newline->is_spanning_u(), newline->const_par_));
		else if(newline->touches(element_[i]))
			element_[i]->addPartialLine(newline);
	}
	}
	#if 0
	this->generateIDs();
	int element_size = element_.size(); 
	for(uint i=0; i<element_size; i++) {
		if(newline->splits(element_[i])) {
			for(int j=0; j<element_[i]->nBasisFunctions(); j++) {
				Basisfunction *f = element_[i]->supportFunction(j);
				if(!newline->containedIn(f) && newline->splits(f)) {
					nNewFunctions += split( const_u, element_[i]->supportFunction(j)->getId(), const_par );
					nRemovedFunctions++;
					j--; // split destroys a Basisfunction which updates the elementpointers
					     // and subsequently removes one basisfunction from the support list
					this->generateIDs();
				}
			}
			element_.push_back(element_[i]->split(!newline->is_spanning_u(), newline->const_par_));
		} else if(newline->touches(element_[i])) {
			element_[i]->addPartialLine(newline);
		}
	}
	#endif
	}

	{ // STEP 2: test every NEW function against ALL old meshlines
#ifdef TIME_LRSPLINE
	PROFILE("STEP 2");
#endif
	meshline_.push_back(newline);
	std::vector<Meshline*>::iterator mit;
	for(uint i=nOldFunctions-nRemovedFunctions-1; i<basis_.size(); i++) {
		for(mit=basis_[i]->partialLineBegin(); mit!=basis_[i]->partialLineEnd(); mit++) {
			if((*mit)->splits(basis_[i]) && 
			   !(*mit)->containedIn(basis_[i])) {
				split( !(*mit)->is_spanning_u(), i, (*mit)->const_par_ );
				i--; // splitting deletes a basisfunction in the middle of the basis_ vector
				break;
			}
		}
		/*
		for(uint j=0; j<meshline_.size(); j++) {
			if(meshline_[j]->splits(basis_[i]) && 
			   !meshline_[j]->containedIn(basis_[i])) {
				split( !meshline_[j]->is_spanning_u(), i, meshline_[j]->const_par_ );
				i--; // splitting deletes a basisfunction in the middle of the basis_ vector
			}
		}
		*/
	}
	}
}

void LRSplineSurface::insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity) {
	insert_line(false, v, start_u, stop_u, multiplicity);
}

int LRSplineSurface::split(bool insert_in_u, int function_index, double new_knot) {
#ifdef TIME_LRSPLINE
	PROFILE("split()");
#endif

	// create the new functions b1 and b2
	Basisfunction b = *basis_[function_index];
	Basisfunction *b1, *b2;
	double *knot = (insert_in_u) ? b.knot_u_  : b.knot_v_;
	int     p    = (insert_in_u) ? b.order_u_ : b.order_v_;
	int     insert_index = 0;
	if(new_knot < knot[0] || knot[p] < new_knot)
		return 0;
	while(knot[insert_index] < new_knot)
		insert_index++;
	double alpha1 = (insert_index == p)  ? 1.0 : (new_knot-knot[0])/(knot[p-1]-knot[0]);
	double alpha2 = (insert_index == 1 ) ? 1.0 : (knot[p]-new_knot)/(knot[p]-knot[1]);
	double newKnot[p+2];
	std::copy(knot, knot+p+1, newKnot+1);
	newKnot[0] = new_knot;
	std::sort(newKnot, newKnot + p+2);
	if(insert_in_u) {
		b1 = new Basisfunction(newKnot  , b.knot_v_, b.controlpoint_, b.dim_, b.order_u_, b.order_v_, b.weight_*alpha1);
		b2 = new Basisfunction(newKnot+1, b.knot_v_, b.controlpoint_, b.dim_, b.order_u_, b.order_v_, b.weight_*alpha2);
	} else { // insert in v
		b1 = new Basisfunction(b.knot_u_, newKnot  , b.controlpoint_, b.dim_, b.order_u_, b.order_v_, b.weight_*alpha1);
		b2 = new Basisfunction(b.knot_u_, newKnot+1, b.controlpoint_, b.dim_, b.order_u_, b.order_v_, b.weight_*alpha2);
	}

	// search for existing function (to make search local b1 and b2 is contained in the *el element)
	std::vector<Element*>::iterator el_it;
	Element *el=NULL;
	for(el_it=b.supportedElementBegin(); el_it<b.supportedElementEnd(); el_it++) {
		if(b1->overlaps(*el_it) && b2->overlaps(*el_it)) {
			el = *el_it;
			break;
		}
	}
	std::vector<Basisfunction*>::iterator it;
	for(it=el->supportBegin(); it != el->supportEnd(); it++) {
		// existing functions need only update their control points and weights
		if(b1 && *b1 == **it) {
			**it += *b1;
			delete b1;
			b1 = NULL;
		} else if(b2 && *b2 == **it) {
			**it += *b2;
			delete b2;
			b2 = NULL;
		}
	}

	// add any brand new functions and detect their support elements
	int newFunctions = 0;
	if(b1) {
		newFunctions++;
		basis_.push_back(b1);
		updateSupport(b1, b.supportedElementBegin(), b.supportedElementEnd());
		b1->inheritEdgeTag(&b, !insert_in_u, true);
		b1->inheritPartialLine(&b);
	}
	if(b2) {
		newFunctions++;
		basis_.push_back(b2);
		updateSupport(b2, b.supportedElementBegin(), b.supportedElementEnd());
		b2->inheritEdgeTag(&b, !insert_in_u, false);
		b2->inheritPartialLine(&b);
	}
	basis_.erase(basis_.begin() + function_index);
	return newFunctions;
}

void LRSplineSurface::getEdgeFunctions(std::vector<Basisfunction*> &edgeFunctions, parameterEdge edge, bool corner) const {
	edgeFunctions.clear();
	for(uint i=0; i<basis_.size(); i++)
		if( (!corner && (edge &  basis_[i]->getEdgeIndex()))  ||
		    ( corner && (edge == basis_[i]->getEdgeIndex()))  ) 
			edgeFunctions.push_back(basis_[i]);
}

void LRSplineSurface::getGlobalKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
	getGlobalUniqueKnotVector(knot_u, knot_v);

	// add in duplicates where apropriate
	for(uint i=0; i<knot_u.size(); i++) {
		for(uint j=0; j<meshline_.size(); j++) {
			if(!meshline_[j]->is_spanning_u() && meshline_[j]->const_par_==knot_u[i]) {
				for(int m=1; m<meshline_[j]->multiplicity_; m++) {
					knot_u.insert(knot_u.begin()+i, knot_u[i]);
					i++;
				}
				break;
			}
		}
	}
	for(uint i=0; i<knot_v.size(); i++) {
		for(uint j=0; j<meshline_.size(); j++) {
			if(meshline_[j]->is_spanning_u() && meshline_[j]->const_par_==knot_v[i]) {
				for(int m=1; m<meshline_[j]->multiplicity_; m++) {
					knot_v.insert(knot_v.begin()+i, knot_v[i]);
					i++;
				}
				break;
			}
		}
	}
}

void LRSplineSurface::getGlobalUniqueKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
	knot_u.clear();
	knot_v.clear();
	// create a huge list of all line instances
	for(uint i=0; i<meshline_.size(); i++) {
		if(meshline_[i]->is_spanning_u())
			knot_v.push_back(meshline_[i]->const_par_);
		else
			knot_u.push_back(meshline_[i]->const_par_);
	}

	// sort them and remove duplicates
	std::vector<double>::iterator it;
	std::sort(knot_u.begin(), knot_u.end());
	std::sort(knot_v.begin(), knot_v.end());
	it = std::unique(knot_u.begin(), knot_u.end());
	knot_u.resize(it-knot_u.begin());
	it = std::unique(knot_v.begin(), knot_v.end());
	knot_v.resize(it-knot_v.begin());
}

void LRSplineSurface::updateSupport(Basisfunction *f,
	                                std::vector<Element*>::iterator start,
	                                std::vector<Element*>::iterator end ) {
#ifdef TIME_LRSPLINE
	PROFILE("update support");
#endif
	std::vector<Element*>::iterator it;
	for(it=start; it!=end; it++)
		if(f->addSupport(*it)) // this tests for overlapping as well as updating  
			(*it)->addSupportFunction(f);
}

void LRSplineSurface::updateSupport(Basisfunction *f) {
	updateSupport(f, element_.begin(), element_.end());
}

void LRSplineSurface::generateIDs() const {
	for(uint i=0; i<basis_.size(); i++) 
		basis_[i]->setId(i);
	for(uint i=0; i<element_.size(); i++) 
		element_[i]->setId(i);
}

bool LRSplineSurface::isLinearIndepByMappingMatrix(bool verbose) const {
#ifdef TIME_LRSPLINE
	PROFILE("Linear independent)");
#endif
	// try and figure out this thing by the projection matrix C

	std::vector<double> knots_u, knots_v;
	getGlobalKnotVector(knots_u, knots_v);
	int nmb_bas = basis_.size();
	int n1 = knots_u.size() - order_u_;
	int n2 = knots_v.size() - order_v_;
	int fullDim = n1*n2;
	bool fullVerbose   = fullDim < 30  && nmb_bas < 50;
	bool sparseVerbose = fullDim < 250 && nmb_bas < 100;

	std::vector<std::vector<boost::rational<long long> > > C;  // rational projection matrix 
	std::vector<double> locKnotU, locKnotV;               // local INDEX knot vectors

	// scaling factor to ensure that all knots are integers (assuming all multiplum of smallest knot span)
	double smallKnotU = 1e300;
	double smallKnotV = 1e300;
	for(uint i=0; i<knots_u.size()-1; i++)
		if(knots_u[i+1]-knots_u[i] < smallKnotU && knots_u[i+1] != knots_u[i])
			smallKnotU = knots_u[i+1]-knots_u[i];
	for(uint i=0; i<knots_v.size()-1; i++)
		if(knots_v[i+1]-knots_v[i] < smallKnotV && knots_v[i+1] != knots_v[i])
			smallKnotV = knots_v[i+1]-knots_v[i];
	double eps = (smallKnotU<smallKnotV) ? smallKnotU/1000 : smallKnotV/1000;

	for (int i = 0; i < nmb_bas; ++i) {
		int startU, startV;
		std::vector<double> locKnotU(basis_[i]->knot_u_, basis_[i]->knot_u_ + basis_[i]->order_u_+1);
		std::vector<double> locKnotV(basis_[i]->knot_v_, basis_[i]->knot_v_ + basis_[i]->order_v_+1);
		
		for(startU=knots_u.size(); startU-->0; )
			if(knots_u[startU] == basis_[i]->knot_u_[0]) break;
		for(int j=0; j<basis_[i]->order_u_; j++) {
			if(knots_u[startU] == basis_[i]->knot_u_[j]) startU--;
			else break;
		}
		startU++;
		for(startV=knots_v.size(); startV-->0; )
			if(knots_v[startV] == basis_[i]->knot_v_[0]) break;
		for(int j=0; j<basis_[i]->order_v_; j++) {
			if(knots_v[startV] == basis_[i]->knot_v_[j]) startV--;
			else break;
		}
		startV++;

		std::vector<boost::rational<long long> > rowU(1,1), rowV(1,1);
		int curU = startU+1;
		for(uint j=0; j<locKnotU.size()-1; j++, curU++) {
			if(locKnotU[j+1] != knots_u[curU]) {
				std::vector<boost::rational<long long> > newRowU(rowU.size()+1, boost::rational<long long>(0));
				for(uint k=0; k<rowU.size(); k++) {
					#define U(x) ((int) (locKnotU[x+k]/smallKnotU + eps))
					int z = (int) (knots_u[curU] / smallKnotU + eps);
					int p = order_u_-1;
					if(z < U(0) || z > U(p+1)) {
						newRowU[k] = boost::rational<long long>(1);
						continue;
					}
					boost::rational<long long> alpha1 = (U(p) <=  z  ) ? 1 : boost::rational<long long>(   z    - U(0),  U(p)  - U(0));
					boost::rational<long long> alpha2 = (z    <= U(1)) ? 1 : boost::rational<long long>( U(p+1) - z   , U(p+1) - U(1));
					newRowU[k]   += rowU[k]*alpha1;
					newRowU[k+1] += rowU[k]*alpha2;
					#undef U
				}
				locKnotU.insert(locKnotU.begin()+(j+1), knots_u[curU]);
				rowU= newRowU;
			}
		}
		int curV = startV+1;
		for(uint j=0; j<locKnotV.size()-1; j++, curV++) {
			if(locKnotV[j+1] != knots_v[curV]) {
				std::vector<boost::rational<long long> > newRowV(rowV.size()+1, boost::rational<long long>(0));
				for(uint k=0; k<rowV.size(); k++) {
					#define V(x) ((int) (locKnotV[x+k]/smallKnotV + eps))
					int z = (int) (knots_v[curV] / smallKnotV + eps);
					int p = order_v_-1;
					if(z < V(0) || z > V(p+1)) {
						newRowV[k] = boost::rational<long long>(1);
						continue;
					}
					boost::rational<long long> alpha1 = (V(p) <=  z  ) ? 1 : boost::rational<long long>(   z    - V(0),  V(p)  - V(0));
					boost::rational<long long> alpha2 = (z    <= V(1)) ? 1 : boost::rational<long long>( V(p+1) - z   , V(p+1) - V(1));
					newRowV[k]   += rowV[k]*alpha1;
					newRowV[k+1] += rowV[k]*alpha2;
					#undef V
				}
				locKnotV.insert(locKnotV.begin()+(j+1), knots_v[curV]);
				rowV= newRowV;
			}
		}
		std::vector<boost::rational<long long> > totalRow(fullDim, boost::rational<long long>(0));
		for(uint i1=0; i1<rowU.size(); i1++)
			for(uint i2=0; i2<rowV.size(); i2++)
				totalRow[(startV+i2)*n1 + (startU+i1)] = rowV[i2]*rowU[i1];

		C.push_back(totalRow);
	}

	if(verbose && sparseVerbose) {
		for(uint i=0; i<C.size(); i++) {
			std::cout << "|";
			for(uint j=0; j<C[i].size(); j++)
				if(C[i][j]==0) {
					if(fullVerbose)
						std::cout << "\t";
					else if(sparseVerbose)
						std::cout << " ";
				} else {
					if(fullVerbose)
						std::cout << C[i][j] << "\t";
					else if(sparseVerbose)
						std::cout << "x";
				}
			std::cout << "|\n";
		}
		std::cout << std::endl;
	}

	// gauss-jordan elimination
	int zeroColumns = 0;
	for(uint i=0; i<C.size() && i+zeroColumns<C[i].size(); i++) {
		boost::rational<long long> maxPivot = 0;
		int maxI = -1;
		for(uint j=i; j<C.size(); j++) {
			if(abs(C[j][i+zeroColumns]) > maxPivot) {
				maxPivot = abs(C[j][i+zeroColumns]);
				maxI = j;
			}
		}
		if(maxI == -1) {
			i--;
			zeroColumns++;
			continue;
		}
		std::vector<boost::rational<long long> > tmp = C[i];
		C[i] = C[maxI];
		C[maxI] = tmp;
		// printf("swapping row %d and %d\n", i, maxI);
		// for(int j=i+zeroColumns+1; j<C[i].size(); j++)
			// C[i][j] /= C[i][i+zeroColumns];
		// C[i][i+zeroColumns] = 1.0;
		for(uint j=i+1; j<C.size(); j++) {
			// if(j==i) continue;
			boost::rational<long long> scale =  C[j][i+zeroColumns] / C[i][i+zeroColumns];
			if(scale != 0) {
				// std::cout << scale << " x row(" << i << ") added to row " << j << std::endl;
				for(uint k=i+zeroColumns; k<C[j].size(); k++)
					C[j][k] -= C[i][k] * scale;
			}
		}
	}

	if(verbose && sparseVerbose) {
		for(uint i=0; i<C.size(); i++) {
			std::cout << "|";
			for(uint j=0; j<C[i].size(); j++)
				if(C[i][j]==0) {
					if(fullVerbose)
						std::cout << "\t";
					else if(sparseVerbose)
				 	std::cout << " ";
				} else {
					if(fullVerbose)
						std::cout << C[i][j] << "\t";
					else if(sparseVerbose)
						std::cout << "x";
				}
			std::cout << "|\n";
		}
		std::cout << std::endl;
	}

	int rank = (nmb_bas < n1*n2-zeroColumns) ? nmb_bas : n1*n2-zeroColumns;
	if(verbose) {
		std::cout << "Matrix size : " << nmb_bas << " x " << n1*n2 << std::endl;
		std::cout << "Matrix rank : " << rank << std::endl;
	}
	
	return rank == nmb_bas;
}

void LRSplineSurface::read(std::istream &is) {
	// first get rid of comments and spaces
	ws(is); 
	char firstChar;
	char buffer[1024];
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	// read actual parameters
	int nBasis, nElements, nMeshlines;
	is >> order_u_;    ws(is);
	is >> order_v_;    ws(is);
	is >> nBasis;      ws(is);
	is >> nMeshlines;  ws(is);
	is >> nElements;   ws(is);
	is >> dim_;        ws(is);
	is >> rational_;   ws(is);
	
	basis_.resize(nBasis);
	meshline_.resize(nMeshlines);
	element_.resize(nElements);

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	// read all basisfunctions
	for(int i=0; i<nBasis; i++) {
		basis_[i] = new Basisfunction(dim_, order_u_, order_v_);
		basis_[i]->read(is);
	}

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	for(int i=0; i<nMeshlines; i++) {
		meshline_[i] = new Meshline();
		meshline_[i]->read(is);
	}

	// get rid of more comments and spaces
	firstChar = is.peek();
	while(firstChar == '#') {
		is.getline(buffer, 1024);
		ws(is);
		firstChar = is.peek();
	}

	for(int i=0; i<nElements; i++) {
		element_[i] = new Element();
		element_[i]->read(is);
		element_[i]->updateBasisPointers(basis_);
	}

}

void LRSplineSurface::write(std::ostream &os) const {
	generateIDs();
	os << "# LRSPLINE\n";
	os << "#\tp1\tp2\tNbasis\tNline\tNel\tdim\trat\n\t";
	os << order_u_ << "\t";
	os << order_v_ << "\t";
	os << basis_.size() << "\t";
	os << meshline_.size() << "\t";
	os << element_.size() << "\t";
	os << dim_ << "\t";
	os << rational_ << "\n";

	std::vector<Element*>::const_iterator eit;
	std::vector<Basisfunction*>::const_iterator bit;
	std::vector<Meshline*>::const_iterator mit;
	os << "# Basis functions:\n";
	int i=0; 
	for(bit=basis_.begin(); bit<basis_.end(); bit++, i++)
		os << **bit << std::endl;
	i = 0;
	os << "# Mesh lines:\n";
	for(mit=meshline_.begin(); mit<meshline_.end(); mit++, i++)
		os << **mit << std::endl;
	os << "# Elements:\n";
	for(eit=element_.begin(); eit<element_.end(); eit++, i++)
		os << **eit << std::endl;
}

void LRSplineSurface::writePostscriptMesh(std::ostream &out) const {
#ifdef TIME_LRSPLINE
	PROFILE("Write EPS");
#endif
	std::vector<double> knot_u, knot_v;
	getGlobalUniqueKnotVector(knot_u, knot_v);
	double min_span_u = knot_u[1] - knot_u[0];
	double min_span_v = knot_v[1] - knot_v[0];
	for(uint i=1; i<knot_u.size()-1; i++)
		min_span_u = (min_span_u<knot_u[i+1]-knot_u[i]) ? min_span_u : knot_u[i+1]-knot_u[i];
	for(uint i=1; i<knot_v.size()-1; i++)
		min_span_v = (min_span_v<knot_v[i+1]-knot_v[i]) ? min_span_v : knot_v[i+1]-knot_v[i];

	// get date
	time_t t = time(0);
	tm* lt = localtime(&t);
	char date[11];
	sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

	// get bounding box
	int dx = end_u_ - start_u_;
	int dy = end_v_ - start_v_;
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;
	// set the duplicate-knot-line (dkl) display width
	double dkl_range = (min_span_u>min_span_v) ? min_span_v*scale/6.0 : min_span_u*scale/6.0; 
	int xmin = (start_u_ - dx/100.0)*scale;
	int ymin = (start_v_ - dy/100.0)*scale;
	int xmax = (end_u_   + dx/100.0)*scale + dkl_range;
	int ymax = (end_v_   + dy/100.0)*scale + dkl_range;

	// print eps header
	out << "%!PS-Adobe-3.0 EPSF-3.0\n";
	out << "%%Creator: LRSplineHelpers.cpp object\n";
	out << "%%Title: LR-spline index domain\n";
	out << "%%CreationDate: " << date << std::endl;
	out << "%%Origin: 0 0\n";
	out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;

/****   USED FOR COLORING OF THE DIAGONAL ELEMENTS SUBJECT TO REFINEMENT (DIAGONAL TEST CASE) ******
	std::vector<BoundingBox> elements = lr->getSegments();
	Point high, low;
	out << ".5 setgray\n";
	for(uint i=0; i<colored.size(); i++) {
		low  = elements[colored[i]].low();
		high = elements[colored[i]].high();

		out << "newpath\n";
		out <<  draw_i_pos[f->i]*scale           << " " <<  draw_j_pos[f->j]*scale << " moveto\n";
		out << draw_i_pos[(f->i+f->width)]*scale << " " <<  draw_j_pos[f->j]*scale << " lineto\n";
		out << draw_i_pos[(f->i+f->width)]*scale << " " << draw_j_pos[(f->j+f->height)]*scale << " lineto\n";
		out <<  draw_i_pos[f->i]*scale           << " " << draw_j_pos[(f->j+f->height)]*scale << " lineto\n";
		out << "closepath\n";
		out << "fill\n";
	}
******************************************************************************************************/

	out << "0 setgray\n";
	out << "1 setlinewidth\n";
	for(uint i=0; i<meshline_.size(); i++) {
		out << "newpath\n";
		double dm = (meshline_[i]->multiplicity_==1) ? 0 : dkl_range/(meshline_[i]->multiplicity_-1);
		for(int m=0; m<meshline_[i]->multiplicity_; m++) {
			if(meshline_[i]->is_spanning_u()) {
				out << meshline_[i]->start_*scale << " " << meshline_[i]->const_par_*scale + dm*m << " moveto\n";
				if(meshline_[i]->stop_ == end_u_)
					out << meshline_[i]->stop_*scale+dkl_range << " " << meshline_[i]->const_par_*scale + dm*m << " lineto\n";
				else
					out << meshline_[i]->stop_*scale << " " << meshline_[i]->const_par_*scale + dm*m << " lineto\n";
			} else {
				out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->start_*scale << " moveto\n";
				if(meshline_[i]->stop_ == end_v_)
					out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->stop_ *scale+dkl_range << " lineto\n";
				else
					out << meshline_[i]->const_par_*scale + dm*m << " " << meshline_[i]->stop_ *scale << " lineto\n";
			}
		}
		out << "stroke\n";
	}

	out << "%%EOF\n";
}

void LRSplineSurface::printElements(std::ostream &out) const {
	for(uint i=0; i<element_.size(); i++) {
		if(i<10) out << " ";
		out << i << ": " << *element_[i] << std::endl;
	}
}

} // end namespace LR

