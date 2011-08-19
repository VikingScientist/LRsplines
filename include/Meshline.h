#ifndef MESHLINE_H
#define MESHLINE_H

namespace LR {

class Meshline {

public:
	Meshline(bool u_line, double const_par, double start, double stop, int multiplicity);

private:
	bool u_line_;
	double const_par_;
	double start_;
	double stop_;
	int multiplicity_;

};

} // end namespace LR

#endif

