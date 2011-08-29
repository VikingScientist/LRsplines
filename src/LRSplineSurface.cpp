
#include "LRSplineSurface.h"
#include "Basisfunction.h"
#include "Meshline.h"
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

typedef unsigned int uint;

namespace LR {


LRSplineSurface::LRSplineSurface(int n1, int n2, int order_u, int order_v, double *knot_u, double *knot_v, double *coef, int dim, bool rational) {

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
			if(i < order_u)   basis_[k-1]->addEdge(WEST);
			if(i>=n1-order_u) basis_[k-1]->addEdge(EAST);
			if(j < order_v)   basis_[k-1]->addEdge(SOUTH);
			if(j>=n2-order_v) basis_[k-1]->addEdge(NORTH);
		}
	}
	for(int i=0; i<n1+order_u; i++) {// const u, spanning v
		int mult = 1;
		while(i<n1+order_u && knot_u[i]==knot_u[i+1]) {
			i++;
			mult++;
		}
		meshline_.push_back(new Meshline(false, knot_u[i], knot_v[0], knot_v[n2], mult) );
	}
	for(int i=0; i<n2+order_v; i++) {// const v, spanning u
		int mult = 1;
		while(i<n2+order_v && knot_v[i]==knot_v[i+1]) {
			i++;
			mult++;
		}
		meshline_.push_back(new Meshline(true, knot_v[i], knot_u[0], knot_u[n2], mult) );
	}
	// element_;
}

void LRSplineSurface::point(Go::Point &pt, double u, double v) const {
	Go::Point cp;
	double basis_ev;
	pt.resize(dim_);
	pt *= 0;
	for(uint i=0; i<basis_.size(); i++) {
		basis_[i]->getControlPoint(cp);
		basis_ev = basis_[i]->evaluate(u,v);
		pt += basis_ev*cp;
	}
}

void LRSplineSurface::point(std::vector<Go::Point> &pts, double u, double v, int derivs) const {
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
		basis_[i]->evaluate(basis_ev, u,v, derivs);
		for(uint j=0; j<pts.size(); j++)
			pts[j] += basis_ev[j]*cp;
	}
}

void LRSplineSurface::computeBasis (double param_u, double param_v, Go::BasisDerivsSf & result ) const {
	result.prepareDerivs(param_u, param_v, 0, -1, basis_.size());
	std::vector<double> values;
	for(uint i=0; i<basis_.size(); i++) {
		basis_[i]->evaluate(values, param_u, param_v, 1);
		result.basisValues[i]   = values[0];
		result.basisDerivs_u[i] = values[1];
		result.basisDerivs_v[i] = values[2];
	}
}


void LRSplineSurface::computeBasis(double param_u, double param_v, Go::BasisPtsSf & result ) const {
	result.preparePts(param_u, param_v, 0, -1, basis_.size());
	for(uint i=0; i<basis_.size(); i++)
		result.basisValues[i] = basis_[i]->evaluate(param_u, param_v);
}

void LRSplineSurface::insert_const_u_edge(double u, double start_v, double stop_v, int multiplicity) {
	insert_line(true, u, start_v, stop_v, multiplicity);
}

void LRSplineSurface::insert_line(bool const_u, double const_par, double start, double stop, int multiplicity) {
	Meshline *newline = new Meshline(!const_u, const_par, start, stop, multiplicity);
	if(multiplicity != 1) {
		std::cerr << "ERROR: LRSplineSurface::insert_const_u_edge() not supported for multiplicity != 1\n";
		std::cerr << "       requested line: " << *newline << std::endl;
	}
	for(uint i=0; i<meshline_.size(); i++) {
		if( meshline_[i]->is_spanning_u() != const_u && meshline_[i]->const_par_ == const_par && 
		     meshline_[i]->start_ <= stop && meshline_[i]->stop_ >= start)  {
			std::cerr << "Meshline extension not implemented yet - just new lines currently\n";
			std::cerr << "New line requeste:    " << *newline      << std::endl;
			std::cerr << "Old line intersected: " << *meshline_[i] << std::endl;
		}
	}
	int nOldFunctions     = basis_.size();  // number of functions before any new ones are inserted
	int nRemovedFunctions = 0;
	int nNewFunctions     = 0; // number of newly created functions eligable for testing in STEP 2

	// STEP 1: test EVERY function against the NEW meshline
	for(int i=0; i<nOldFunctions-nRemovedFunctions; i++) {
		if(newline->splits(basis_[i])) {
			std::cout << "Meshline " << *newline << " splits function " << *basis_[i] << std::endl;
			nNewFunctions += split( const_u, i, const_par );
			i--; // splitting deletes a basisfunction in the middle of the basis_ vector
			nRemovedFunctions++;
		}
	}

	// STEP 2: test every NEW function against ALL old meshlines
	meshline_.push_back(newline);
	for(uint i=nOldFunctions-nRemovedFunctions-1; i<basis_.size(); i++) {
		for(uint j=0; j<meshline_.size(); j++) {
			if(meshline_[j]->splits(basis_[i]) && 
			   !meshline_[j]->containedIn(basis_[i])) {
				std::cout << "Old meshline " << *meshline_[j] << " splits function " << *basis_[i] << std::endl;
				split( !meshline_[j]->is_spanning_u(), i, meshline_[j]->const_par_ );
				i--; // splitting deletes a basisfunction in the middle of the basis_ vector
			}
		}
	}
}

void LRSplineSurface::insert_const_v_edge(double v, double start_u, double stop_u, int multiplicity) {
	insert_line(false, v, start_u, stop_u, multiplicity);
}

int LRSplineSurface::split(bool insert_in_u, int function_index, double new_knot) {
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
	for(uint i=0; i<basis_.size(); i++) {
		if(b1 && *b1 == *basis_[i]) {
			std::cout << "Basis " << i << " equals new basis b1\n";
			std::cout << "(" << *basis_[i] << ") == (" << *b1 << ")\n";
			*basis_[i] += *b1;
			delete b1;
			b1 = NULL;
		} else if(b2 && *b2 == *basis_[i]) {
			std::cout << "Basis " << i << " equals new basis b2\n";
			std::cout << "(" << *basis_[i] << ") == (" << *b2 << ")\n";
			*basis_[i] += *b2;
			delete b2;
			b2 = NULL;
		}
	}
	int newFunctions = 0;
	if(b1) newFunctions++;
	if(b2) newFunctions++;
	if(b1) basis_.push_back(b1);
	if(b2) basis_.push_back(b2);
	basis_.erase(basis_.begin() + function_index);
	return newFunctions;
}

void LRSplineSurface::getGlobalKnotVector(std::vector<double> &knot_u, std::vector<double> &knot_v) const {
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

void LRSplineSurface::read(std::istream &is) {
}

void LRSplineSurface::write(std::ostream &os) const {
	std::vector<Basisfunction*>::const_iterator bit;
	std::vector<Meshline*>::const_iterator mit;
	os << "# Basis functions:\n";
	int i=0; 
	for(bit=basis_.begin(); bit<basis_.end(); bit++, i++)
		os << i << ": " << **bit << std::endl;
	i = 0;
	os << "# Mesh lines:\n";
	for(mit=meshline_.begin(); mit<meshline_.end(); mit++, i++)
		os << i << ": " << **mit << std::endl;
}

void LRSplineSurface::writePostscriptMesh(std::ostream &out) {
	std::vector<double> knot_u, knot_v;
	getGlobalKnotVector(knot_u, knot_v);
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

} // end namespace LR

