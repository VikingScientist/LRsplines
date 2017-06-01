#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Profiler.h"
#include "LRSpline/Element.h"
#include "LRSpline/Meshline.h"

using namespace LR;
using namespace std;

vector<double> bounding_box(const LRSplineSurface &surf) {
	vector<double> result(4);
	result[0] =  1e9;
	result[1] = -1e9;
	result[2] =  1e9;
	result[3] = -1e9;
	for(Basisfunction *b : surf.getAllBasisfunctions()) {
		std::vector<double>::const_iterator cp = b->cp();
		result[0] = (cp[0] < result[0]) ? cp[0] : result[0];
		result[1] = (cp[0] > result[1]) ? cp[0] : result[1];
		result[2] = (cp[1] < result[2]) ? cp[1] : result[2];
		result[3] = (cp[1] > result[3]) ? cp[1] : result[3];
	}
	return result;
}

void writePostscriptElements(std::ostream &out, const LRSplineSurface &lr, int nu, int nv, double scale, double color[3]) {
	// Set color
	out << color[0] << " ";
	out << color[1] << " ";
	out << color[2] << " ";
	out << "setrgbcolor \n";
	for(int iEl=0; iEl<lr.nElements(); iEl++) {
		double umin = lr.getElement(iEl)->umin();
		double umax = lr.getElement(iEl)->umax();
		double vmin = lr.getElement(iEl)->vmin();
		double vmax = lr.getElement(iEl)->vmax();

		std::vector<double> pt;
		lr.point(pt, umin, vmin, iEl);
		out << "newpath\n";
		out <<  pt[0]*scale << " " << pt[1]*scale << " moveto\n";
		for(int i=1; i<nu; i++) { // SOUTH side
			double u = umin + (umax-umin)*i/(nu-1);
			double v = vmin;
			lr.point(pt, u, v, iEl, false, true);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=1; i<nv; i++) { // EAST side
			double u = umax;
			double v = vmin + (vmax-vmin)*i/(nv-1);
			lr.point(pt, u, v, iEl, false, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=nu-1; i-->0; ) { // NORTH side
			double u = umin + (umax-umin)*i/(nu-1);
			double v = vmax;
			lr.point(pt, u, v, iEl, true, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=nv-1; i-->1; ) { // WEST side
			double u = umin;
			double v = vmin + (vmax-vmin)*i/(nv-1);
			lr.point(pt, u, v, iEl, true, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		out << "closepath\n";
		out << "fill\n";
	}

	out << "0 setgray\n";
	out << "1 setlinewidth\n";

	for(int iEl=0; iEl<lr.nElements(); iEl++) {
		double umin = lr.getElement(iEl)->umin();
		double umax = lr.getElement(iEl)->umax();
		double vmin = lr.getElement(iEl)->vmin();
		double vmax = lr.getElement(iEl)->vmax();

		std::vector<double> pt;
		lr.point(pt, umin, vmin, iEl);
		out << "newpath\n";
		out <<  pt[0]*scale << " " << pt[1]*scale << " moveto\n";
		for(int i=1; i<nu; i++) { // SOUTH side
			double u = umin + (umax-umin)*i/(nu-1);
			double v = vmin;
			lr.point(pt, u, v, iEl, false, true);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=1; i<nv; i++) { // EAST side
			double u = umax;
			double v = vmin + (vmax-vmin)*i/(nv-1);
			lr.point(pt, u, v, iEl, false, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=nu-1; i-->0; ) { // NORTH side
			double u = umin + (umax-umin)*i/(nu-1);
			double v = vmax;
			lr.point(pt, u, v, iEl, true, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		for(int i=nv-1; i-->1; ) { // WEST side
			double u = umin;
			double v = vmin + (vmax-vmin)*i/(nv-1);
			lr.point(pt, u, v, iEl, true, false);
			out <<  pt[0]*scale << " " << pt[1]*scale << " lineto\n";
		}
		out << "closepath\n";
		out << "stroke\n";
	}

}

int main(int argc, char **argv) {
	// error test input
	if(argc < 2) {
		cout << "Usage: " << argv[0] << " <inputfiles>\n";
		exit(1);
	}

	// set default patch colors
//int color[][3] = {{  0  ,114  ,190},
//                  {218  , 83  , 25},
//                  {238  ,178  , 32},
//                  {126  , 47  ,142},
//                  {119  ,173  , 48},
//                  { 77  ,191  ,239},
//                  {163  , 20  , 47} };
	double color[][3] = {{     0 ,  0.4470 ,  0.7410},
	                     {0.8500 ,  0.3250 ,  0.0980},
	                     {0.9290 ,  0.6940 ,  0.1250},
	                     {0.4940 ,  0.1840 ,  0.5560},
	                     {0.4660 ,  0.6740 ,  0.1880},
	                     {0.3010 ,  0.7450 ,  0.9330},
	                     {0.6350 ,  0.0780 ,  0.1840}};


	// read all patches
	vector<LRSplineSurface*> surfs;
	for(int i=1; i<argc; ++i) {
		ifstream inputfile;
		inputfile.open(argv[i]);
		if(!inputfile.is_open()) {
			cerr << "Error: could not open file " << argv[i] << endl;
			exit(2);
		}
		while(inputfile.is_open() && !inputfile.eof()) {
			LRSplineSurface *lr = new LRSplineSurface();
			inputfile >> *lr;
			surfs.push_back(lr);
		}
	}

	// compute global bounding box
	vector<double> bb;
	double x[2];
	double y[2];
	x[0] =  1e9;
	x[1] = -1e9;
	y[0] =  1e9;
	y[1] = -1e9;
	for(auto srf : surfs) {
		bb = bounding_box(*srf);
		x[0] = min(bb[0], x[0]);
		x[1] = max(bb[1], x[1]);
		y[0] = min(bb[2], y[0]);
		y[1] = max(bb[3], y[1]);
	}

	// get date
	time_t t = time(0);
	tm* lt = localtime(&t);
	char date[11];
	sprintf(date, "%02d/%02d/%04d", lt->tm_mday, lt->tm_mon + 1, lt->tm_year+1900);

	// compute scaling factors
	double dx = x[1]-x[0];
	double dy = y[1]-y[0];
	double scale = (dx>dy) ? 1000.0/dx : 1000.0/dy;
	int xmin = (x[0] - dx/20.0)*scale;
	int ymin = (y[0] - dy/20.0)*scale;
	int xmax = (x[1] + dx/20.0)*scale;
	int ymax = (y[1] + dy/20.0)*scale;

	// open output file
	ofstream out("mesh.eps");

	// print eps header
	out << "%!PS-Adobe-3.0 EPSF-3.0\n";
	out << "%%Creator: LRSplineHelpers.cpp object\n";
	out << "%%Title: LR-spline physical domain\n";
	out << "%%CreationDate: " << date << std::endl;
	out << "%%Origin: 0 0\n";
	out << "%%BoundingBox: " << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl;

	// draw all patches
	int i=0;
	for(auto srf : surfs) {
		writePostscriptElements(out, *srf, 7,7, scale, color[i++ % 7]);
	}

	// close file
	out << "%%EOF\n";

	out.close();

  cout << "Written mesh.eps" << endl;

}
